library(data.table)
library(caret)
library(BGLR)
library(asreml)
library(ASRgenomics)


set.seed(12345)
sink("../../bayescaccuracy/yuanshilog.txt")
M012 = fread("../../raw/all.raw")

G = as.matrix(M012[,!c(1:6)])
rownames(G) = M012$FID

data <- asreml.read.table("../../survivaladjust.csv",sep=",",header=T)
gid = data$ID
G <- G[match(gid, rownames(G)), ]
G[1:10,1:10]
dim(G)

# 清洗基因组数据：maf，callrate
G_clean = qc.filtering(G,marker.callrate = 0.05,ind.callrate = 0.1)
clean_m012 = G_clean$M.clean
dim(clean_m012)
sum(is.na(clean_m012))


clean_m012bayes <- clean_m012[rownames(clean_m012) %in% gid, ]
dim(clean_m012bayes)

X = clean_m012

# 封装为函数
run_cv_analysis <- function(data, seed_values = 1:10) {
  # 筛选非NA数据
  data1 <- data[!is.na(data$survival), ]
  phedur <- data1
  
  bayesre_final <- list()
  
  for(iii in seed_values) {
    set.seed(iii)
    w <- createFolds(1:length(phedur$ID), k = 5)
    
    # 创建5个折叠数据集
    fold_data <- list()
    for(i in 1:5) {
      d <- phedur
      d[w[[i]], 4] <- NA  # 将第i折的数据设为NA
      fold_data[[i]] <- d
    }
    
    # 用于存储每个折叠的结果
    fold_results <- list()
    
    # 对每个折叠进行GBLUP分析
    for(fold_idx in 1:5) {
      d <- fold_data[[fold_idx]]
      
      # 设置数据类型和因子
      setDF(d)
      for(i in 1:3) d[, i] <- as.factor(d[, i])
      
      phedurdbayesa = BGLR(d$survival,ETA = list(list(~d$Tank*d$Sex,
                                          X=clean_m012bayes,model="BayesC")),
                    nIter = 5000,burnIn = 500,thin = 5,saveAt = "",
                    S0 = NULL,R2 = 0.5, weights = NULL,verbose = TRUE,
                    rmExistingFiles = TRUE, groups=NULL)
      # 提取随机效应
      snpd = phedurdbayesa$ETA[[1]]$b
      ghat = X %*% snpd
      bad = data.frame(ID = row.names(ghat), effect = ghat[,1])
      # 获取预测值和真实值
      na_rows <- d[!complete.cases(d), ]
      badselect <- merge(bad,na_rows,by="ID")
      badselect <- badselect[,c(1,2)]
      badblup <- merge(badselect,phedur,by="ID")
      # 计算准确性和偏差
      accurbayesa = cor(badblup$adjust,badblup$effect) # 准确性
      unbiabayesa = coef(lm(badblup$adjust~badblup$effect))[[2]]
      fold_results[[fold_idx]] <- data.frame(Accur = accurbayesa, Unbia = unbiabayesa)            
      }
    
    # 合并5个折叠的结果
    refinal <- do.call(rbind, fold_results)
    
    # 计算标准误的函数
    se <- function(x) {
      sd(x) / sqrt(length(x))
    }
    
    # 计算平均值和标准误
    accur_mean <- mean(refinal$Accur)
    accur_se <- se(refinal$Accur)
    unbia_mean <- mean(refinal$Unbia)
    unbia_se <- se(refinal$Unbia)
    
    # 存储该seed的结果
    bayesre_final[[iii]] <- data.frame(
      accur = accur_mean,
      accur_se = accur_se,
      unbia = unbia_mean,
      unbia_se = unbia_se
    )
  }
  # 计算所有seed的平均结果
  accur_values <- sapply(bayesre_final, function(x) x$accur)
  accur_se_values <- sapply(bayesre_final, function(x) x$accur_se)
  unbia_values <- sapply(bayesre_final, function(x) x$unbia)
  unbia_se_values <- sapply(bayesre_final, function(x) x$unbia_se)
  
  mean_accur <- mean(accur_values)
  mean_accur_se <- mean(accur_se_values)
  mean_unbia <- mean(unbia_values)
  mean_unbia_se <- mean(unbia_se_values)
  
  # 返回最终结果
  result <- data.frame(
    accur = mean_accur,
    accur_se = mean_accur_se,
    unbia = mean_unbia,
    unbia_se = mean_unbia_se
  )
  print(bayesre_final)
  return(result)
}

# 设置循环次数
n_loops <- 1
results_list <- list()

for(loop in 1:n_loops) {
  cat(paste("\n=== 第", loop, "次循环 ===\n"))
  
  # 复制原始数据
  data_missing <- data
  
  # 设置不同的随机种子，确保每次缺失的位置不同
  set.seed(1000 + loop * 100)  # 使用不同的种子
  
  # 随机选择1000个索引设为缺失
  missing_indices <- sample(1:nrow(data), 0)
  
  # 将这些行的survival设为NA
  data_missing[missing_indices, "survival"] <- NA
  
   
  # 检查缺失情况
  missing_count <- sum(is.na(data$survival))
  cat(paste("设置了", missing_count, "个缺失值\n"))
  cat(paste("缺失比例:", round(missing_count/nrow(data)*100, 2), "%\n"))
  
  # 运行分析
  result <- run_cv_analysis(data_missing)
  result$loop <- loop
  result$missing_count <- missing_count
  result$missing_indices_seed <- 1000 + loop * 100
  
  # 保存结果
  results_list[[loop]] <- result

}

# 显示结果

print("平均值:")
print(results_list)

sink()
q()



