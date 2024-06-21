library(org.At.tair.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

data <- read.table("test_3\\TCGA_HNSC_DGE_GENE.txt")
data = data[order(data[,'padj'], decreasing = TRUE),]
head(data)
table(data$sig)
up_geneList <- bitr(rownames(data[data$sig == 'up',]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
down_geneList <- bitr(rownames(data[data$sig == 'down',]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
diff_geneList <- bitr(rownames(data), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(up_geneList)
head(down_geneList)
head(diff_geneList)
kegg_diff = enrichKEGG(gene = diff_geneList$ENTREZID, organism = 'hsa',keyType = 'kegg',
                       pvalueCutoff = 0.99,qvalueCutoff = 0.99)
kegg_up = enrichKEGG(gene = up_geneList$ENTREZID, organism = 'hsa',keyType = 'kegg',
                     pvalueCutoff = 0.99, qvalueCutoff = 0.99)
kegg_down = enrichKEGG(gene = down_geneList$ENTREZID, organism = 'hsa',keyType = 'kegg',
                       pvalueCutoff = 0.99, qvalueCutoff = 0.99)
search_kegg_organism('hsa',by = 'kegg_code')
#此处先省去id转换
down_kegg = kegg_down[kegg_down$p.adjust < 0.05,]
down_kegg$group = -1

up_kegg = kegg_up[kegg_up$p.adjust < 0.05,]
up_kegg$group = 1

dat = rbind(up_kegg, down_kegg)
dat$p.adjust = log10(dat$p.adjust)
dat$p.adjust = -dat$p.adjust
dat = dat[order(-dat$p.adjust),]
dat$p.adjust = dat$p.adjust*dat$group
dat$group = ifelse(dat$group == -1, 'down', 'up')
dat$group = factor(as.character(dat$group),levels = c('down','up'))

top10_dat = dat[1:10,]
top10_dat = top10_dat[order(top10_dat$p.adjust),]
g_kegg = ggplot(top10_dat, aes(x = factor(Description,levels = Description),
                         y = p.adjust, fill = group)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#37DBD0','#FF7C6C')) +
  scale_x_discrete(name = 'Top 10 Pathway names') +
  scale_y_continuous(name = '(±)|log10P-value|') + coord_flip() +
  ggtitle('Top 10 Pathway Enrichment')
ggsave("test_3\\g_Top10kegg_bar.png", g_kegg, width = 10, height = 8,dpi = 300)
ggsave('test_3\\g_Top10kegg_bar.pdf', g_kegg, width = 10, height = 8)

#ggplot2可视化：按显著性和富集因子排序的气泡图
library(stringr)
top10_dat$GeneRatio = sapply(str_split(top10_dat$GeneRatio, '/'),function(x){
  as.numeric(x)[1] / as.numeric(x)[2]
})
g_dot <-
  ggplot(top10_dat,aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = Count, fill = p.adjust),
             alpha = 1,
             shape = 21,
             stroke=0.5,
  ) + 
  scale_x_continuous(expand = c(0.05,0)) +
  scale_fill_gradient(low="#1062c9",high="#e1415f",name="±|p|") +
  labs(color=expression(p.adjust),
       size="Count",
       x="GeneRatio",
       y=" ") +
  guides(fill = NULL,
         size = guide_legend(order=1)) +
  theme_linedraw(base_size = 14) +
  theme(axis.text = element_text(color="black"),
        axis.ticks.length = unit(0.25,"cm")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", linewidth = 0.25),
        panel.border = element_rect(fill=NA,linewidth = 1),
        legend.background = element_blank()
  )
g_dot
ggsave("test_3\\g_Top10kegg_dot.png", g_dot, width = 10, height = 8,dpi = 300)
ggsave('test_3\\g_Top10kegg_dot.pdf', g_dot, width = 10, height = 8)
