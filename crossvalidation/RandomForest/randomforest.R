library(randomForest)
library(data.table)
library(ASRgenomics)
library(asreml)
library(caret)


sink("../RandomForestaccuracy/allSNP.txt")

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

# 使用之前找到的超参数
grid_results <- read.csv("grid_results_RF.csv",header = T,sep = ",")
best_by_rmse <- grid_results[which.min(grid_results$mean_rmse), ]

print(best_by_rmse)

#进行交叉验证
n_repeats <- 10
n_folds <- 5

# 存储所有结果
all_repeats <- list()

for(rep in 1:n_repeats) {
  cat("          第", rep, "次重复\n")
  
  set.seed(1000 + rep)  # 每次重复用不同种子
  cv_folds <- createFolds(1:length(y), k = n_folds, list = TRUE)
  
  # 存储本次重复的5折结果
  repeat_results <- data.frame()
  
  for(f in 1:n_folds) {
    # 分割训练集和测试集
    test_idx <- cv_folds[[f]]
    train_idx <- setdiff(1:length(y), test_idx)
    
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    

    # 训练模型（使用固定超参数）
    set.seed(789 + f)
    
    # 处理maxnodes参数
    maxnodes_val <- if(is.null(best_by_rmse$maxnodes) || best_by_rmse$maxnodes == "NULL") NULL else as.numeric(best_by_rmse$maxnodes)
    
    rf_model <- randomForest(
      x = X_train,
      y = y_train,
      ntree = best_by_rmse$ntree,
      mtry = best_by_rmse$mtry,
      nodesize = best_by_rmse$nodesize,
      maxnodes = maxnodes_val,
      importance = FALSE,
      keep.forest = TRUE
    )
    
    # 预测
    y_pred <- predict(rf_model, X_test)
    
    # 计算指标
    cor_val <- cor(y_test, y_pred)
    rmse_val <- sqrt(mean((y_test - y_pred)^2))
    bia_val <- coef(lm(y_test~y_pred))[[2]]
    
    
    # 记录本折结果
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
  repeat_summary <- data.frame(
    repetition = rep,
    mean_cor = mean(repeat_results$correlation),
    sd_cor = sd(repeat_results$correlation),
    mean_rmse = mean(repeat_results$rmse),
    sd_rmse = sd(repeat_results$rmse),
    mean_bias = mean(repeat_results$bia),
    sd_bias = sd(repeat_results$bia)
    
  )
  
  all_repeats[[rep]] <- list(
    fold_results = repeat_results,
    summary = repeat_summary
  )
  
  cat("\n第", rep, "次重复平均 - 相关系数:", 
      round(repeat_summary$mean_cor, 4), 
      "±", round(repeat_summary$sd_cor, 4), "\n")
}


# 汇总所有重复的结果

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("          10次重复5折交叉验证结果汇总\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 提取所有重复的汇总结果
all_summaries <- do.call(rbind, lapply(all_repeats, function(x) x$summary))
print(all_summaries)


# 计算总体统计
final_summary <- data.frame(
  metric = c("correlation", "rmse", "bia"),
  mean = c(
    mean(all_summaries$mean_cor),
    mean(all_summaries$mean_rmse),
    mean(all_summaries$mean_bias)
  ),
  se = c(
    sd(all_summaries$mean_cor)/sqrt(nrow(all_summaries)),
    sd(all_summaries$mean_rmse)/sqrt(nrow(all_summaries)),
    sd(all_summaries$mean_bias)/sqrt(nrow(all_summaries))
  )
)

print(final_summary)
cat("结束时间", format(Sys.time(), "%H:%M:%S"), "\n")
sink()