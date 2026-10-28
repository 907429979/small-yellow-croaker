library(randomForest)
library(data.table)
library(ASRgenomics)
library(asreml)
library(caret)
library(parallel)
library(doParallel)
library(foreach)

sink("searchHyperparameter.txt")
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("           随机森林 基因组选择 - 搜索超参数\n")
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

# 5. 训练随机森林模型 ========================================================

# 定义需要网格搜索的参数
param_grid <- expand.grid(
  ntree = c(300, 500, 800, 1000),                 # 树的数量
  mtry = c(floor(sqrt(ncol(X_train_val)) * 0.5),   # 特征采样数 - 50%
           floor(sqrt(ncol(X_train_val))),          # 特征采样数 - 100% (默认)
           floor(sqrt(ncol(X_train_val)) * 1.5),    # 特征采样数 - 150%
           floor(ncol(X_train_val) * 0.1)),         # 特征采样数 - 10% 的SNP
  nodesize = c(5, 10, 15, 20),                      # 叶节点最小样本数
  maxnodes = c(50, 100, NULL)                       # 最大叶节点数（NULL表示无限制）
)
cat("参数组合总数:", nrow(param_grid), "\n")

set.seed(456)
n_search <- min(100, nrow(param_grid))  # 随机森林训练较慢，适当减少搜索数量
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

# 使用foreach进行并行计算
grid_results <- foreach(i = 1:nrow(param_grid), .combine = rbind, 
                        .packages = c("randomForest", "caret")) %dopar% {
                          
                          params <- param_grid[i, ]
                          
                          fold_rmse <- c()
                          fold_cor <- c()
                          
                          for(f in 1:n_folds) {
                            # 分割训练验证集内的训练和验证
                            inner_train_idx <- setdiff(1:length(y_train_val), folds[[f]])
                            inner_valid_idx <- folds[[f]]
                            
                            X_inner_train <- X_train_val[inner_train_idx, ]
                            y_inner_train <- y_train_val[inner_train_idx]
                            X_inner_valid <- X_train_val[inner_valid_idx, ]
                            y_inner_valid <- y_train_val[inner_valid_idx]
                            
                            # 训练随机森林模型
                            set.seed(789 + f)
                            
                            # 处理maxnodes参数
                            maxnodes_val <- if(is.null(params$maxnodes) || params$maxnodes == "NULL") NULL else as.numeric(params$maxnodes)
                            
                            rf_model <- randomForest(
                              x = X_inner_train,
                              y = y_inner_train,
                              ntree = params$ntree,
                              mtry = params$mtry,
                              nodesize = params$nodesize,
                              maxnodes = maxnodes_val,
                              importance = FALSE,
                              keep.forest = TRUE
                            )
                            
                            # 预测
                            pred <- predict(rf_model, X_inner_valid)
                            
                            # 计算指标
                            rmse <- sqrt(mean((y_inner_valid - pred)^2))
                            cor_val <- ifelse(sd(pred) > 0, cor(y_inner_valid, pred), 0)
                            
                            fold_rmse <- c(fold_rmse, rmse)
                            fold_cor <- c(fold_cor, cor_val)
                          }
                          
                          # 返回结果
                          data.frame(
                            ntree = params$ntree,
                            mtry = params$mtry,
                            nodesize = params$nodesize,
                            maxnodes = if(is.null(params$maxnodes)) "NULL" else as.character(params$maxnodes),
                            mean_rmse = mean(fold_rmse),
                            sd_rmse = sd(fold_rmse),
                            mean_cor = mean(fold_cor),
                            sd_cor = sd(fold_cor),
                            stringsAsFactors = FALSE
                          )
                        }

stopCluster(cl)

# 保存网格搜索结果
write.csv(grid_results, file = "grid_results_RF.csv", row.names = FALSE)

# 找到最优参数 ======================================
grid_results
# 基于RMSE
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]

cat("\n搜索超参数结束时间", format(Sys.time(), "%H:%M:%S"), "\n")

# 选择最优参数用于后续建模 
best_params <- best_by_rmse
cat("\n=== 选择的最优参数 ===\n")
print(best_params)

cat("\n=== 用最优参数训练最终模型 ===\n")

# 处理maxnodes参数
maxnodes_final <- if(best_params$maxnodes == "NULL") NULL else as.numeric(best_params$maxnodes)

# 训练最终模型
set.seed(101)
final_model <- randomForest(
  x = X_train_val,
  y = y_train_val,
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  maxnodes = maxnodes_final,
  importance = FALSE,
  keep.forest = TRUE
)

# 在测试集上评估
y_pred <- predict(final_model, X_test)
test_cor <- cor(y_test, y_pred)
test_rmse <- sqrt(mean((y_test - y_pred)^2))
test_bias <- coef(lm(y_test ~ y_pred))[2]

cat("\n=== 测试集评估结果 ===\n")
cat("测试集相关系数:", round(test_cor, 4), "\n")
cat("测试集RMSE:", round(test_rmse, 4), "\n")
cat("测试集偏差:", round(test_bias, 4), "\n")

# 保存测试集数据
ceshijidata <- survival_data[-train_idx, ]
write.csv(ceshijidata, file = "ceshijidata_RF.csv", row.names = FALSE)


sink()
