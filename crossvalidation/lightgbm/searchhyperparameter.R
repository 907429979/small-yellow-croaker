library(lightgbm)
library(data.table)
library(ASRgenomics)
library(caret)
library(parallel)
library(doParallel)
library(foreach)
sink("searchHyperparameter.txt")
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("           LightGBM 基因组选择 - 搜索超参数\n")
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
data <- read.csv("../survivaladjust.csv", sep = ",", header = TRUE)

for(i in 1:3) data[, i] <- as.factor(data[, i])


survival_data <- data[!is.na(data$bw), ]
cat("bw有效观测数:", nrow(survival_data), "\n")



# 4. 构建特征矩阵 ============================================================

X <- clean_m012[match(survival_data$ID, rownames(clean_m012)), ]
y <- survival_data$adjust

set.seed(123)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train_val <- X[train_idx, ]  # 训练+验证集（用于调优）
y_train_val <- y[train_idx]
X_test <- X[-train_idx, ]      # 独立测试集（用于最终评估）
y_test <- y[-train_idx]

cat("训练验证集大小:", length(y_train_val), "\n")
cat("独立测试集大小:", length(y_test), "\n")
cat("测试集比例:", round(length(y_test)/length(y)*100, 2), "%\n")

# 5. 训练LightGBM模型 ========================================================

# 创建LightGBM数据集
dtrain <- lgb.Dataset(X, label = y)

#定义需要网格搜索的参数

param_grid <- expand.grid(
  num_leaves = c(31, 63, 95, 127),           # 树复杂度
  learning_rate = c(0.01, 0.03, 0.05, 0.1),  # 学习率
  min_data_in_leaf = c(20, 50, 100),          # 叶子节点最小样本数
  feature_fraction = c(0.3, 0.5, 0.7),        # SNP采样比例
  bagging_fraction = c(0.6, 0.7, 0.8, 0.9),   # 每棵树使用的样本比例
  bagging_freq = c(1, 3, 5, 10),              # bagging的频率
  lambda_l1 = c(0, 0.1, 1),                    # L1正则化
  lambda_l2 = c(0, 0.1, 1)                     # L2正则化
)
cat("参数组合总数:", nrow(param_grid), "\n")

set.seed(456)
n_search <- min(100, nrow(param_grid))
search_idx <- sample(1:nrow(param_grid), n_search)
param_grid <- param_grid[search_idx, ]
cat("实际搜索组合数:", nrow(param_grid), "\n")


# 6. 在训练验证集上做交叉验证调优 ============================================

set.seed(789)
n_folds <- 5
# 只在训练验证集上创建交叉验证折
folds <- createFolds(1:length(y_train_val), k = n_folds, list = TRUE)

# 设置并行计算
n_cores <- detectCores() - 64
cl <- makeCluster(n_cores)

registerDoParallel(cl)
cat("使用", n_cores, "核并行计算\n")

# 存储结果
grid_results <- data.frame()

grid_results <- foreach(i = 1:nrow(param_grid), 
                        .combine = rbind,
                        .packages = c("lightgbm", "caret")) %dopar% {
                          
                          params <- param_grid[i, ]
                          
                          fold_rmse <- c()
                          fold_cor <- c()
                          
                          for(f in 1:n_folds) {
                            # 分割数据
                            inner_train_idx <- setdiff(1:length(y_train_val), folds[[f]])
                            inner_valid_idx <- folds[[f]]
                            
                            X_inner_train <- X_train_val[inner_train_idx, ]
                            y_inner_train <- y_train_val[inner_train_idx]
                            X_inner_valid <- X_train_val[inner_valid_idx, ]
                            y_inner_valid <- y_train_val[inner_valid_idx]
                            
                            # 创建数据集
                            dtrain <- lgb.Dataset(X_inner_train, label = y_inner_train)
                            
                            # 设置参数（每个核心内用单线程）
                            lgb_params <- list(
                              objective = "regression",
                              metric = "rmse",
                              boosting_type = "gbdt",
                              num_leaves = as.integer(params$num_leaves),
                              learning_rate = params$learning_rate,
                              feature_fraction = params$feature_fraction,
                              bagging_fraction = params$bagging_fraction,
                              bagging_freq = params$bagging_freq,
                              min_data_in_leaf = as.integer(params$min_data_in_leaf),
                              lambda_l1 = params$lambda_l1,
                              lambda_l2 = params$lambda_l2,
                              verbose = -1,
                              num_threads = 1,
                              seed = 789 + f
                            )
                            
                            # 训练模型
                            model <- lgb.train(
                              params = lgb_params,
                              data = dtrain,
                              nrounds = 500,
                              verbose = 0
                            )
                            
                            # 预测和评估
                            pred <- predict(model, X_inner_valid)
                            rmse <- sqrt(mean((y_inner_valid - pred)^2))
                            cor_val <- ifelse(sd(pred) > 0, cor(y_inner_valid, pred), 0)
                            
                            fold_rmse <- c(fold_rmse, rmse)
                            fold_cor <- c(fold_cor, cor_val)
                          }
                          
                          # 返回结果
                          data.frame(
                            num_leaves = params$num_leaves,
                            learning_rate = params$learning_rate,
                            feature_fraction = params$feature_fraction,
                            bagging_fraction = params$bagging_fraction,
                            bagging_freq = params$bagging_freq,
                            min_data_in_leaf = params$min_data_in_leaf,
                            lambda_l1 = params$lambda_l1,
                            lambda_l2 = params$lambda_l2,
                            mean_rmse = mean(fold_rmse),
                            sd_rmse = sd(fold_rmse),
                            mean_cor = mean(fold_cor),
                            sd_cor = sd(fold_cor)
                          )
                        }

stopCluster(cl)

write.csv(grid_results,file="grid_results.csv",row.names = FALSE)

#找到最优参数======================================

# 基于RMSE
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]
print(best_by_rmse)
cat("搜素超参数结束时间", format(Sys.time(), "%H:%M:%S"), "\n")

#选择最优参数用于后续建模 
best_params <- best_by_rmse
cat("\n=== 选择的最优参数 ===\n")
print(best_params)

cat("\n=== 用最优参数训练最终模型 ===\n")

dtrain_full <- lgb.Dataset(X_train_val, label = y_train_val)

final_params <- list(
  objective = "regression",
  metric = "rmse",
  boosting_type = "gbdt",
  num_leaves = as.integer(best_params$num_leaves),
  learning_rate = best_params$learning_rate,
  feature_fraction = best_params$feature_fraction,
  bagging_fraction = best_params$bagging_fraction,
  bagging_freq = as.integer(best_params$bagging_freq),
  min_data_in_leaf = best_params$min_data_in_leaf,
  lambda_l1 = best_params$lambda_l1,
  lambda_l2 = best_params$lambda_l2,
  verbose = 1,
  num_threads = n_cores,
  seed = 101
)

final_model <- lgb.train(
  params = final_params,
  data = dtrain_full,
  nrounds = 500,
  verbose = 1
)

y_pred <- predict(final_model, X_test)
test_cor <- cor(y_test, y_pred)
test_rmse <- sqrt(mean((y_test - y_pred)^2))
test_bias <- coef(lm(y_test~y_pred))[[2]]

cat("测试集相关系数:", round(test_cor, 4), "\n")
cat("测试集RMSE:", round(test_rmse, 4), "\n")
ceshijidata <- survival_data[!data$ID %in% train_idx,]
write.csv(ceshijidata,file = "ceshijidata.csv",row.names = FALSE)

sink()
