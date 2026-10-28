#GO富集
library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)
library(clusterProfiler)

zhushi <- fread("../genetic/all_gene_function.xls")
gogene <- zhushi[,c(1,9)]                          #将基因ID和GOID单独成列
gogenefilter <- subset(gogene, GOs != "_")
gotermlong <- gogenefilter %>%
  separate_rows(GOs, sep = ";") %>%
  # 去除可能的空格
  mutate(GOs = trimws(GOs)) %>%
  # 过滤掉空值（如果有）
  filter(GOs != "")

#搜索所有的GO描述
go_id <- gotermlong$GOs
go_id <- unique(go_id)

godescription <- data.frame(
  go_id = go_id
)

#寻找GOterm
library(GO.db)

godescription <- select(GO.db, 
              keys = godescription$go_id,      # 你的GO ID向量
              columns = c("TERM", "ONTOLOGY"),       # 我们想查询TERM（描述）
              keytype = "GOID")       # 指定keys的类型是GOID

write.csv(godescription,file="../genetic/godescription.csv",row.names = FALSE)



colnames(godescription)[1] <-"go_id"
colnames(gotermlong) <- c("gene_id","go_id")

bggo <- left_join(gotermlong,godescription,by="go_id")
bggo <- na.omit(bggo)

write.csv(bggo,file="../enrichanalyze/goenrich/bggo.csv",row.names = FALSE)

geneid <- gene$geneID
colnames(bggo)[3] <- "Description"
enrich_result2 <- enricher(
  gene      = geneid,
  TERM2GENE = bggo[,c(2,1)],
  TERM2NAME = bggo[,c(2,3)],
  pvalueCutoff  = 0.05,
  pAdjustMethod = "holm",
  qvalueCutoff  = 0.2
)


head(enrich_result)

enrich_result@result$ONTOLOGY <- bggo$ONTOLOGY[
  match(enrich_result@result$ID, bggo$go_id)
]

result2 <- enrich_result2@result
write.csv(result,file="./goenrich/allenrichresult.csv",row.names = FALSE)



# 可视化：条形图


keggfilter <- keggresult[keggresult$pvalue < 0.05,]
keggfilter$logp <- -log10(keggfilter$pvalue)

# 按 Count 排序（让条形从上到下递增）
kegg_df <- keggfilter %>% arrange(Count)


ggplot(kegg_df, aes(x = Count, y = reorder(Description, Count), fill = logp)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradientn(
    colours = c("blue", "purple", "red"),  # 蓝→紫→红
    name = expression(-log[10](P.value))   # 下标10
  ) +
  labs(x = "Gene count", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(family = "serif",size = 14, colour = "black"),
    axis.text.x = element_text(family = "serif",size = 14, colour = "black"),
    axis.title.x = element_text(family = "serif",size = 14),
    legend.title = element_text(family = "serif",size = 14),
    legend.text = element_text(family = "serif",size = 14),
    panel.grid.major.y = element_blank(),  # 移除横向网格线
    panel.grid.minor = element_blank()
  )




barplot(statuenrich_result, showCategory = 20, title = "GO Enrichment Analysis")

# 可视化：气泡图
goresult <- result[result$pvalue < 0.05,]

# 计算 -log10(pvalue)
goresult$logp <- -log10(goresult$pvalue)

# 筛选每类前N个
library(dplyr)
top_n <- 10
go10filter <- goresult %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = pvalue, n = top_n) %>%
  ungroup()

# ggplot2 手动绘制
ggplot(go10filter, aes(x = Count, y = reorder(Description, Count))) +
  geom_point(aes(size = Count, color = logp)) +
  scale_color_gradientn(
    colours = c("#313695",  # 深蓝
                "#4575B4",  # 蓝
                "#74ADD1",  # 浅蓝
                "#FDAE61",  # 橙
                "#F46D43",  # 红橙
                "#D73027",  # 红
                "#A50026"), # 深红,
    name = "-log10(pvalue)"
  ) +
  scale_size_continuous(range = c(2, 6), name = "count") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Gene Count", y = NULL) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))



ggsave("keggenrichment.tiff", 
       plot = last_plot(),    # 或指定ggplot对象，如 plot = p
       device = "tiff",
       dpi = 600,
       width = 10, 
       height = 8, 
       units = "in",
       compression = "lzw")   # 压缩减少文件大小



