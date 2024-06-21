#1、加载包
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

###----
#基因id转换，kegg和go富集用的ID类型是ENTREZID）
data <- read.table("test_3\\TCGA_HNSC_DGE_GENE.txt")
head(data)

diff <- bitr(rownames(data), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- diff$ENTREZID

#GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#4、将结果保存到当前路径
ego_results_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)

write.csv(ego_ALL,file = "test_3\\ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "test_3\\ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "test_3\\ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "test_3\\ego_result_MF.csv",row.names = T)

##将以上通路重新组合成带因子的数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  GeneRatio=c(ego_result_BP$GeneRatio, ego_result_CC$GeneRatio, ego_result_MF$GeneRatio),
  p.adjust=c(ego_result_BP$p.adjust, ego_result_CC$p.adjust, ego_result_MF$p.adjust),
  type=factor(c(rep("biological process", length(row.names(ego_result_BP))), 
                rep("cellular component", length(row.names(ego_result_CC))),
                rep("molecular function", length(row.names(ego_result_MF)))), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
### ----
# 直接导入本地数据验证


### top10 GO term和top10通路
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ego = go_enrich_df
#Count
ego = ego[order(ego$GeneNumber, decreasing = TRUE),]
head(ego$GeneNumber)
max(ego$GeneNumber)
factor(ego[1:10,]$Description, levels = ego[1:10,]$Description)
top10count_GO_term_bar = ggplot(data=ego[1:10,], 
                                aes(x=factor(ego[1:10,]$Description, levels = ego[1:10,]$Description),y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Top 10 Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
top10count_GO_term_bar
ggsave("test_3\\top10GeneNumber_GO_term_bar.png", top10count_GO_term_bar, width = 10, height = 8,dpi = 300)
ggsave('test_3\\top10GeneNumber_GO_term_bar.pdf', top10count_GO_term_bar, width = 10, height = 8)


library(stringr)
ego$GeneRatio = sapply(str_split(ego$GeneRatio, '/'),function(x){
  as.numeric(x)[1] / as.numeric(x)[2]
})
top10pvalue_GO_term_dot = ggplot(ego[1:10,],
       aes(y=Description,x=GeneRatio))+
  geom_point(aes(size=GeneNumber,color=-log10(p.adjust)))+
  scale_color_gradient(low = "#f0faff", high = "#1062c9")+
  labs(color=expression(-log10(p.adjust),size="GeneNumber"), 
       x="Gene Ratio",y="GO term",title="Top 10 GO Enrichment")+
  theme_bw()
top10pvalue_GO_term_dot
ggsave("test_3\\top10pvalue_GO_term_dot.png", top10pvalue_GO_term_dot, width = 10, height = 8,dpi = 300)
ggsave('test_3\\top10pvalue_GO_term_dot.pdf', top10pvalue_GO_term_dot, width = 10, height = 8)
