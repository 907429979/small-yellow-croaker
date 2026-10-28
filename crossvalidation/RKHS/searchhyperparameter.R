library(BGLR)
library(data.table)
library(ASRgenomics)
library(asreml)
library(caret)
library(parallel)
library(doParallel)

sink("searchHyperparameter_RKHS.txt")
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("BGLR RKHS 基因组选择 - 搜索最优带宽\n")
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

# 构建特征矩阵
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

# ==================== RKHS 带宽参数搜索 ====================

# 计算距离矩阵用于带宽参数选择
cat("\n=== 计算遗传距离矩阵 ===\n")
# 使用标准化后的SNP计算欧氏距离
X_scaled <- scale(X_train_val)
D_sq <- as.matrix(dist(X_scaled))^2
# 获取距离的分位数作为带宽候选的参考
d_quantiles <- quantile(D_sq[D_sq > 0], probs = c(0.1, 0.5, 0.9))
cat("距离矩阵分位数:", round(d_quantiles, 4), "\n")

# 定义带宽参数搜索网格
# 基于距离分位数构建合理的带宽范围
h_grid <- c(
  1 / d_quantiles[3],  # 基于高分位数（较小带宽，局部性更强）
  1 / d_quantiles[2],  # 基于中位数（中等带宽）
  1 / d_quantiles[1],  # 基于低分位数（较大带宽，全局性更强）
  mean(1 / d_quantiles),
  # 添加一些额外的候选值以确保覆盖全面
  0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10
)

# 去重并排序
h_grid <- unique(sort(h_grid[h_grid > 0]))
cat("\n带宽参数搜索网格:\n")
print(round(h_grid, 6))

cat("\n总共有", length(h_grid), "个带宽参数需要搜索\n")

# 设置交叉验证
set.seed(789)
n_folds <- 5
folds <- createFolds(1:length(y_train_val), k = n_folds, list = TRUE)

# 设置并行计算
n_cores <- detectCores() - 64
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("使用", n_cores, "核并行计算\n")

# 存储结果
grid_results <- data.frame()

# 开始带宽参数搜索
for(i in 1:length(h_grid)) {
  h <- h_grid[i]
  
  cat("\n正在搜索带宽 h =", round(h, 6), "(", i, "/", length(h_grid), ")\n")
  
  # 计算当前带宽下的核矩阵（使用训练验证集）
  K <- exp(-h * D_sq)
  
  fold_results <- foreach(f = 1:n_folds, .combine = rbind, 
                          .packages = c("BGLR")) %dopar% {
                            
                            # 分割训练验证集内的训练和验证
                            inner_train_idx <- setdiff(1:length(y_train_val), folds[[f]])
                            inner_valid_idx <- folds[[f]]
                            
                            y_inner_train <- y_train_val[inner_train_idx]
                            y_inner_valid <- y_train_val[inner_valid_idx]

                            
                            # 使用BGLR的预测功能
                            n_train <- length(y_inner_train)
                            n_valid <- length(y_inner_valid)
                            
                            # 构建包含NA的完整响应向量
                            y_all <- c(y_inner_train, rep(NA, n_valid))
                            
                            # 构建完整的核矩阵
                            K_all <- K[c(inner_train_idx, inner_valid_idx), 
                                       c(inner_train_idx, inner_valid_idx)]
                            
                            ETA_all <- list(list(K = K_all, model = "RKHS"))
                            
                            fm_pred <- BGLR(y = y_all,
                                            ETA = ETA_all,
                                            nIter = 2000,
                                            burnIn = 500,
                                            verbose = FALSE,
                                            saveAt = tempfile())
                            
                            pred <- fm_pred$yHat[(n_train+1):(n_train+n_valid)]
                            
                            # 计算指标
                            rmse <- sqrt(mean((y_inner_valid - pred)^2))
                            cor_val <- ifelse(sd(pred) > 0, cor(y_inner_valid, pred), 0)
                            
                            data.frame(rmse = rmse, cor = cor_val)
                          }
  
  # 记录结果
  grid_results <- rbind(grid_results, data.frame(
    bandwidth = h,
    mean_rmse = mean(fold_results$rmse),
    sd_rmse = sd(fold_results$rmse),
    mean_cor = mean(fold_results$cor),
    sd_cor = sd(fold_results$cor)
  ))
  
  # 输出当前进度
  cat("  完成 - 平均RMSE:", round(mean(fold_results$rmse), 4), 
      " 平均相关系数:", round(mean(fold_results$cor), 4), "\n")
  cat("  当前最佳相关系数:", max(grid_results$mean_cor, na.rm = TRUE), 
      " (带宽:", grid_results$bandwidth[which.max(grid_results$mean_cor)], ")\n")
}

stopCluster(cl)

# 保存搜索结果
write.csv(grid_results, file = "RKHS_bandwidth_results.csv", row.names = FALSE)

# 找到最优带宽
best_by_cor <- grid_results[which.max(grid_results$mean_cor), ]
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    带宽搜索完成\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n=== 基于相关系数的最优带宽 ===\n")
print(best_by_cor)

cat("\n=== 基于RMSE的最优带宽 ===\n")
print(best_by_rmse)

cat("结束时间", format(Sys.time(), "%H:%M:%S"), "\n")

# 选择最优带宽（这里选择基于RMSE的）
best_bandwidth <- best_by_rmse$bandwidth
cat("\n=== 选择的最优带宽 ===\n")
cat("最优带宽 (基于RMSE):", round(best_bandwidth, 6), "\n")

# ==================== 用最优带宽训练最终模型 ====================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("          用最优带宽训练最终RKHS模型\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 使用训练验证集训练最终模型
cat("\n计算最终模型的核矩阵...\n")

X_train_scaled <- scale(X_train_val)
X_test_scaled <- scale(X_test, 
                       center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale"))

X_all_scaled <- rbind(X_train_scaled, X_test_scaled)
n_train <- nrow(X_train_scaled)
n_test <- nrow(X_test_scaled)

# 构建完整的核矩阵
K_all <- exp(-best_bandwidth * as.matrix(dist(X_all_scaled))^2)

# 构建包含NA的完整响应向量
y_all <- c(y_train_val, rep(NA, n_test))

ETA_all <- list(list(K = K_all, model = "RKHS"))

# 用完整数据重新训练并预测测试集
fm_test <- BGLR(y = y_all,
                ETA = ETA_all,
                nIter = 2000,
                burnIn = 500,
                verbose = FALSE,
                saveAt = tempfile())

y_pred <- fm_test$yHat[(n_train+1):(n_train+n_test)]

# 评估测试集性能
test_cor <- cor(y_test, y_pred)
test_rmse <- sqrt(mean((y_test - y_pred)^2))
test_bias <- coef(lm(y_test ~ y_pred))[2]

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    测试集评估结果\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("最优带宽:", round(best_bandwidth, 6), "\n")
cat("测试集相关系数:", round(test_cor, 4), "\n")
cat("测试集RMSE:", round(test_rmse, 4), "\n")
cat("测试集回归系数（bias）:", round(test_bias, 4), "\n")
cat("预测值标准差:", round(sd(y_pred), 4), "\n")
cat("观测值标准差:", round(sd(y_test), 4), "\n")

ceshijidata <- survival_data[!data$ID %in% train_idx,]
write.csv(ceshijidata,file = "ceshijidata.csv",row.names = FALSE)

sink()