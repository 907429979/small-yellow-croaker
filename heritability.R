library(asreml)
library(ASRgenomics)
library(data.table)
library(openxlsx)
library(tidyverse)

M012 =  fread("../genetic/all.raw")
M012[1:10,1:10]

G = as.matrix(M012[,!c(1:6)])
rownames(G) = M012$FID
G[1:10,1:10]

# 清洗基因组数据：maf，callrate
G_clean = qc.filtering(G,marker.callrate = 0.1,ind.callrate = 0.1)
clean_m012 = G_clean$M.clean
clean_m012[1:10,1:10]
dim(clean_m012)

Gmat = G.matrix(G_clean$M.clean)$G
Gmat[1:5,1:5]

# G逆矩阵构建，转化为asreml-R需要的的格式
diag(Gmat) = diag(Gmat) + 0.01
ginv = G.inverse(Gmat,sparseform = T)$Ginv.sparse
head(ginv)

##读取表型数据
data <- asreml.read.table("../phe.csv",sep=",",header=T)
data <- as.data.frame(data)
for(i in 1:3) data[,i] = as.factor(data[,i])
str(data)

##构建GBLUP模型
tlmodel = asreml(tl ~ 1+Tank*Sex, 
                     random = ~ vm(ID,ginv), 
                     dense = ~ vm(ID,ginv),
                     residual = ~ idv(units), 
                     workspace = "10Gb",
                     data=data)
wald.asreml(tlmodel)
summary(tlmodel)$varcomp
vpredict(tlmodel, h2 ~ V1/(V1+V2))

blmodel = asreml(bl ~ 1+Tank*Sex, 
                     random = ~ vm(ID,ginv), 
                     dense = ~ vm(ID,ginv),
                     residual = ~ idv(units), 
                     workspace = "10Gb",
                     data=data)
wald.asreml(blmodel)
summary(blmodel)$varcomp
vpredict(blmodel, h2 ~ V1/(V1+V2))

bwmodel = asreml(bw ~ 1+Tank*Sex, 
                     random = ~ vm(ID,ginv), 
                     dense = ~ vm(ID,ginv),
                     residual = ~ idv(units), 
                     workspace = "10Gb",
                     data=data)
wald.asreml(bwmodel)
summary(bwmodel)$varcomp
vpredict(bwmodel, h2 ~ V1/(V1+V2))

phedurmodel = asreml(phe_dur ~ 1+Tank*Sex, 
                         random = ~ vm(ID,ginv), 
                         dense = ~ vm(ID,ginv),
                         residual = ~ idv(units), 
                         workspace = "10Gb",
                         data=data)
wald.asreml(phedurmodel)
summary(phedurmodel)$varcomp
vpredict(phedurmodel, yuanshih2 ~ V1/(V1+V2))

