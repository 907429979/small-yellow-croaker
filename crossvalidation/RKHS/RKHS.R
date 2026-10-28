library(BGLR)
library(data.table)
library(ASRgenomics)
library(asreml)
library(caret)
library(parallel)
library(doParallel)

sink("../RKHSaccuracy/GWAS/6kaccuracy.txt")
cat("          重复10次的五折交叉验证（使用最优带宽）\n")
cat("开始时间", format(Sys.time(), "%H:%M:%S"), "\n")

cat("\n=== 读取基因型数据 ===\n")
M012 <- fread("../raw/GWAS/Rformat/6k.raw")
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

# 使用之前找到的最优带宽
grid_results <- read.csv("grid_results.csv",header = T,sep = ",")
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]

best_bandwidth <- best_by_rmse$bandwidth
cat("最优带宽:", round(best_bandwidth, 6), "\n")

# 设置重复次数和折数
n_repeats <- 10
n_folds <- 5

# 存储所有结果
cv_results_all <- data.frame()

# 开始重复交叉验证
for(rep in 1:n_repeats) {
  cat("\n=== 第", rep, "/", n_repeats, "次重复 ===\n")
  
  # 为每次重复创建新的CV折
  set.seed(1000 + rep)  # 每次重复用不同种子
  folds <- createFolds(1:length(y), k = n_folds, list = TRUE)
  
  # 存储本次重复的结果
  repeat_results <- data.frame()
  
  # 对每一折进行验证
  for(f in 1:n_folds) {
    # 分割训练集和验证集
    train_idx <- setdiff(1:length(y), folds[[f]])
    test_idx <- folds[[f]]
    
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    n_train <- length(y_train)
    n_test <- nrow(X_test)

    # 预测验证集
    # 构建训练集和验证集之间的核矩阵
    X_train_scaled <- scale(X_train)
    X_test_scaled <- scale(X_test, 
                           center = attr(X_train_scaled, "scaled:center"),
                           scale = attr(X_train_scaled, "scaled:scale"))
    X_all_cv <- rbind(X_train_scaled, X_test_scaled)
    D_all_sq <- as.matrix(dist(X_all_cv))^2
    K_all <- exp(-best_bandwidth * D_all_sq)

    
    # 使用BGLR预测验证集
    y_all <- c(y_train, rep(NA, n_test))
    
    ETA_all <- list(list(K = K_all, model = "RKHS"))
    
    fm_pred <- BGLR(y = y_all,
                    ETA = ETA_all,
                    nIter = 5000,
                    burnIn = 500,
                    verbose = FALSE,
                    saveAt = tempfile())
    
    y_pred <- fm_pred$yHat[(n_train+1):(n_train+n_test)]
    
    # 计算验证集指标
    cor_val <- cor(y_test, y_pred)
    rmse_val <- sqrt(mean((y_test - y_pred)^2))
    bia_val <- coef(lm(y_test~y_pred))[[2]]
  
    repeat_results <- rbind(repeat_results, data.frame(
      fold = f,
      correlation = cor_val,
      rmse = rmse_val,
      bia = bia_val,
      n_train = length(train_idx),
      n_test = length(test_idx)
    ))
    
    cat("  第", f, "折 - 相关系数:", round(cor_val, 4), 
        "bias:", round(bia_val, 4), "\n")
  }
  
  # 计算本次重复的平均值
  rep_summary <- data.frame(
    repeat_id = rep,
    mean_rmse = mean(repeat_results$rmse),
    sd_rmse = sd(repeat_results$rmse),
    mean_cor = mean(repeat_results$correlation),
    sd_cor = sd(repeat_results$correlation),
    mean_bias = mean(repeat_results$bia),
    sd_bias = sd(repeat_results$bia)
  )
  
  # 输出本次重复的汇总
  cat("\n第", rep, "次重复汇总:\n")
  cat("  平均RMSE:", round(rep_summary$mean_rmse, 4), 
      "±", round(rep_summary$sd_rmse, 4), "\n")
  cat("  平均相关系数:", round(rep_summary$mean_cor, 4), 
      "±", round(rep_summary$sd_cor, 4), "\n")
  
  # 合并结果
  cv_results_all <- rbind(cv_results_all, rep_summary)
}
cat("          打印10次结果\n")
print(cv_results_all)

# ==================== 交叉验证结果汇总 ====================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("          10×5折交叉验证最终结果\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 计算总体统计


cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("          10次重复5折交叉验证结果汇总\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 计算总体统计
final_summary <- data.frame(
  metric = c("correlation", "rmse", "bia"),
  mean = c(
    mean(cv_results_all$mean_cor),
    mean(cv_results_all$mean_rmse),
    mean(cv_results_all$mean_bias)
  ),
  se = c(
    sd(cv_results_all$mean_cor)/sqrt(nrow(cv_results_all)),
    sd(cv_results_all$mean_rmse)/sqrt(nrow(cv_results_all)),
    sd(cv_results_all$mean_bias)/sqrt(nrow(cv_results_all))
  )
)

print(final_summary)

cat("结束时间", format(Sys.time(), "%H:%M:%S"), "\n")
sink()



