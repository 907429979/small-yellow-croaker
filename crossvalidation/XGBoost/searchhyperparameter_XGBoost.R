library(xgboost)
library(data.table)
library(ASRgenomics)
library(asreml)
library(caret)
library(parallel)
library(doParallel)
library(foreach)

sink("searchHyperparameter_GBDT.txt")
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("           XGBoost 基因组选择 - 搜索超参数\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("开始时间", format(Sys.time(), "%H:%M:%S"), "\n")


cat("\n=== 读取基因型数据 ===\n")
M012 <- fread("../raw/all.raw")
cat("原始基因型数据维度:", dim(M012), "\n")
G <- as.matrix(M012[, !c(1:6)])
rownames(G) <- M012$FID

# 基因型数据质控
G_clean <- qc.filtering(G, marker.callrate = 0.1, ind.callrate = 0.1)
clean_m012 <- G_clean$M.clean
cat("清洗后SNP矩阵维度:", dim(clean_m012), "\n")

# 读取表型数据
cat("\n=== 读取表型数据 ===\n")
data <- asreml.read.table("../survivaladjust.csv", sep = ",", header = TRUE)

for(i in 1:3) data[, i] <- as.factor(data[, i])

survival_data <- data[!is.na(data$survival), ]
cat("survival有效观测数:", nrow(survival_data), "\n")

# ===========================================================================
# 构建特征矩阵
# ===========================================================================
X <- clean_m012[match(survival_data$ID, rownames(clean_m012)), ]
y <- survival_data$adjust

set.seed(123)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train_val <- X[train_idx, ]  # 训练+验证集（用于调优）
y_train_val <- y[train_idx]
X_test <- X[-train_idx, ]      # 独立测试集
y_test <- y[-train_idx]

cat("训练验证集大小:", length(y_train_val), "\n")
cat("独立测试集大小:", length(y_test), "\n")
cat("测试集比例:", round(length(y_test)/length(y)*100, 2), "%\n")

# ===========================================================================
# XGBoost 超参数网格定义
# ===========================================================================
# 针对你的数据特点（945个体，高维SNP）设计的参数范围
param_grid <- expand.grid(
  # 核心参数
  nrounds = c(100, 300, 500, 800),           # 迭代次数（树的数量）
  max_depth = c(3, 4, 5, 6),                  # 树的最大深度（防止过拟合）
  eta = c(0.01, 0.05, 0.1),              # 学习率（越小越稳健）
  
  lambda = c(0, 0.5, 1, 2),                    # L2正则化项
  alpha = c(0, 0.1, 0.5),                       # L1正则化项
  
  # 采样参数
  subsample = c(0.6, 0.8, 1.0),                 # 样本采样比例
  colsample_bytree = c(0.3, 0.5, 0.7),     # 特征采样比例（对高维数据重要）
  
  # 最小叶子节点权重（相当于随机森林的nodesize）
  min_child_weight = c(1, 3, 5, 7)
)

cat("参数组合总数:", nrow(param_grid), "\n")

set.seed(456)
n_search <- min(100, nrow(param_grid))  
search_idx <- sample(1:nrow(param_grid), n_search)
param_grid <- param_grid[search_idx, ]
cat("实际搜索组合数:", nrow(param_grid), "\n")

# ===========================================================================
# 交叉验证设置
# ===========================================================================
set.seed(789)
n_folds <- 5
folds <- createFolds(1:length(y_train_val), k = n_folds, list = TRUE)

# 设置并行计算
n_cores <- detectCores() - 64
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("使用", n_cores, "核并行计算\n")

# ===========================================================================
# 并行搜索超参数
# ===========================================================================

# 存储结果
grid_results <- data.frame()
grid_results <- foreach(i = 1:nrow(param_grid), .combine = rbind,
                        .packages = c("xgboost", "caret")) %dopar% {
                          
                          params <- param_grid[i, ]
                          local_X_train_val <- X_train_val
                          local_y_train_val <- y_train_val
                          local_folds <- folds
                          local_n_folds <- n_folds
                          fold_rmse <- c()
                          fold_cor <- c()
                          
                          for(f in 1:n_folds) {
                              # 分割训练验证集内的训练和验证
                              inner_train_idx <- setdiff(1:length(y_train_val), folds[[f]])
                              inner_valid_idx <- folds[[f]]
                              
                              X_inner_train <- local_X_train_val[inner_train_idx, , drop = FALSE]
                              y_inner_train <- local_y_train_val[inner_train_idx]
                              X_inner_valid <- local_X_train_val[inner_valid_idx, , drop = FALSE]
                              y_inner_valid <- local_y_train_val[inner_valid_idx]
                              
                              # 创建xgb.DMatrix
                              dtrain <- xgb.DMatrix(X_inner_train, label = y_inner_train)
                              dvalid <- xgb.DMatrix(X_inner_valid, label = y_inner_valid)
                              
                              # 设置参数
                              xgb_params <- list(
                                objective = "reg:squarederror",
                                eval_metric = "rmse",
                                booster = "gbtree",
                                max_depth = as.integer(params$max_depth),
                                eta = params$eta,
                                subsample = params$subsample,
                                colsample_bytree = params$colsample_bytree,
                                min_child_weight = params$min_child_weight,
                                lambda = params$lambda,
                                alpha = params$alpha,
                                nthread = 1,                     # 每个线程内不并行，避免冲突
                                seed = 789 + f
                              )
                              
                              # 训练模型（使用watchlist观察验证集）
                              watchlist <- list(train = dtrain, valid = dvalid)
                              model <- xgb.train(
                                params = xgb_params,
                                data = dtrain,
                                nrounds = params$nrounds,
                                watchlist = watchlist,
                                verbose = 0,
                                early_stopping_rounds = 50
                              )
                              
                              # 预测
                              pred <- predict(model, dvalid)
                              
                              # 计算指标
                              rmse <- sqrt(mean((y_inner_valid - pred)^2))
                              cor_val <- ifelse(sd(pred) > 0, cor(y_inner_valid, pred), 0)
                              
                              fold_rmse <- c(fold_rmse, rmse)
                              fold_cor <- c(fold_cor, cor_val)
                            }
                          
                          
                          # 返回结果
                          data.frame(
                            nrounds = params$nrounds,
                            max_depth = params$max_depth,
                            eta = params$eta,
                            subsample = params$subsample,
                            colsample_bytree = params$colsample_bytree,
                            min_child_weight = params$min_child_weight,
                            lambda = params$lambda,
                            alpha = params$alpha,
                            mean_rmse = mean(fold_rmse),
                            sd_rmse = sd(fold_rmse),
                            mean_cor = mean(fold_cor),
                            sd_cor = sd(fold_cor)
                          )
                        }
  
stopCluster(cl)

# ===========================================================================
# 保存和筛选最优参数
# ===========================================================================
write.csv(grid_results, file = "grid_results_XGBoost.csv", row.names = FALSE)

# 基于RMSE找到最优参数
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]

cat("\n搜索超参数结束时间", format(Sys.time(), "%H:%M:%S"), "\n")
cat("结束时间", format(Sys.time(), "%H:%M:%S"), "\n")

cat("\n=== 最优参数（基于RMSE）===\n")
print(best_by_rmse)

# ===========================================================================
# 用最优参数训练最终模型
# ===========================================================================
cat("\n=== 用最优参数训练最终模型 ===\n")

# 准备完整训练集和测试集的DMatrix
dtrain_full <- xgb.DMatrix(data = X_train_val, label = y_train_val)
dtest <- xgb.DMatrix(data = X_test, label = y_test)

# 设置最优参数
final_params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  max_depth = best_by_rmse$max_depth,
  eta = best_by_rmse$eta,
  lambda = best_by_rmse$lambda,
  alpha = best_by_rmse$alpha,
  subsample = best_by_rmse$subsample,
  colsample_bytree = best_by_rmse$colsample_bytree,
  min_child_weight = best_by_rmse$min_child_weight,
  verbosity = 1
)

# 训练最终模型（使用交叉验证找到最优迭代次数）
set.seed(101)
cv_result <- xgb.cv(
  params = final_params,
  data = dtrain_full,
  nrounds = best_by_rmse$nrounds,
  nfold = 5,
  early_stopping_rounds = 50,
  verbose = TRUE,
  prediction = TRUE
)

# 最优迭代次数
best_nrounds_final <- nrow(cv_result$evaluation_log)
cat("\n最优迭代次数:", best_nrounds_final, "\n")

# 用全部训练数据训练最终模型
final_model <- xgb.train(
  params = final_params,
  data = dtrain_full,
  nrounds = best_nrounds_final,
  verbose = 1
)

# ===========================================================================
# 在测试集上评估
# ===========================================================================
y_pred <- predict(final_model, dtest)

test_cor <- cor(y_test, y_pred)
test_rmse <- sqrt(mean((y_test - y_pred)^2))
test_bias <- coef(lm(y_test ~ y_pred))[2]

cat("\n=== 测试集评估结果 ===\n")
cat("测试集相关系数:", round(test_cor, 4), "\n")
cat("测试集RMSE:", round(test_rmse, 4), "\n")
cat("测试集偏差（回归系数）:", round(test_bias, 4), "\n")


ceshijidata <- survival_data[-train_idx, ]
write.csv(ceshijidata, file = "ceshijidata_XGBoost.csv", row.names = FALSE)

cat("结束时间", format(Sys.time(), "%H:%M:%S"), "\n")
sink()
