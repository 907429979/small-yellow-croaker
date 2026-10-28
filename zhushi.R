library(dplyr)
library(readr)
library(data.table)
# 读入显著SNP列表


# 写出为 BED 文件，方便下游 bedtools 等工具使用
sigstatu <- fread("../status/sigstatu.csv")
statuintervals <- sigstatu[,c(2,3,3)]

colnames(statuintervals)[2] <- "Start"
colnames(statuintervals)[3] <- "End"

fwrite(statuintervals, "../status/statu.bed")



test <- statugenezhushi[,c(1,8,9,6)]
result <- sigsurvival %>%
  left_join(test, by = "chrom",relationship = "many-to-many")  %>% filter(pos >= start, pos <= end) 

statugeno <- result[,c(1,4,8,12,10,11)]


 
zhushistatu <- fread("./结果/status/gene.txt")  


gene_rowsstatu <- zhushistatu[zhushistatu$V6 == "gene", ] 
gene_rowsstatu <- gene_rowsstatu[,-c(1,2,3,9,10,11)]
colnames(gene_rowsstatu) <- c("chrom","way","gene","genestart","geneend","geneID")



colnames(test) <- c("chrom","Start","End","ID")



all.gene.functiion <- fread("D:/Asreml/大黄鱼/annotation/all.gene.function.xls")


gene_rowsstatu$geneID <- gsub("ID=", "", gene_rowsstatu$geneID)
gene_rowsstatu$geneID <- gsub(";", "", gene_rowsstatu$geneID)

statugenezhushi <- gene_rowsstatu %>%
  left_join(all.gene.functiion, by = "geneID")
fwrite(statugenezhushi,file="./结果/status/statusgenezhushi.csv")


library(stringr)
result2 <- str_remove(text, "^[^ ]+ ")
result2 <- str_remove(result2, " \\[.+\\]$")

statugenezhushi$name <- str_remove(statugenezhushi$name, " \\[.+\\]$")
statusgenename <- statugenezhushi[,c(1,2,3,4,5,6,15)]
fwrite(statusgenename,file="./结果/status/statugenename.csv")
