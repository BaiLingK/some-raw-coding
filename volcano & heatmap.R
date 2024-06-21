load("test_3/TCGA_HNSC_EXP.Rdata")
exp_gene = gene_counts
head(exp_gene)
# 根据样本名barcode分组>>>
group_list=ifelse(as.numeric(substr(colnames(exp_gene),14,15)) < 10,'tumor','normal')
# 使用DESeq2对肿瘤样本和对照样本进行差异表达分析>>>
library(DESeq2)
coldata <- data.frame(
  sampleID = colnames(exp_gene),
  groups = group_list
)
coldata$groups <- factor(coldata$groups, levels = c('tumor', 'normal'))

round_counts <- round(exp_gene)
dds = DESeqDataSetFromMatrix(
  countData = round_counts,
  colData = coldata,
  design = ~groups)



dds_counts <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = TRUE)
res <- results(dds_counts, contrast = c("groups", "tumor", "normal"))
rownames(res) <- rownames(exp_gene)
head(res)
write.csv(res, file = "test_3\\TCGA_HNSC_DGE.csv")
# 找到显著上调或下调的基因>>>
pvalue_threshold <- 0.01
log2fc_threshold <- 1

res = res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)),]

res[which(res$padj < pvalue_threshold & res$log2FoldChange >= log2fc_threshold), 'sig'] <- 'up'
res[which(res$padj < pvalue_threshold & res$log2FoldChange <= -log2fc_threshold), 'sig'] <- 'down'
res[which((res$padj >= pvalue_threshold) | (abs(res$log2FoldChange) <= log2fc_threshold)), 'sig'] <- 'none'
table(res$sig)
DEG_GENE <- subset(res, sig %in% c('up', 'down'))
write.table(DEG_GENE, file = "test_3\\TCGA_HNSC_DGE_GENE.txt", sep = "\t")

# 绘制火山图>>>
library(ggplot2)

max(res$log2FoldChange)
min(res$log2FoldChange)
max(-log10(res$padj))
min(-log10(res$padj))

p = ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +
  scale_color_manual(values = c('#ff5f41', 'gray77','#54ff1d'), limits = c('up','none','down')) +
  labs(x = 'log2 Fold Change', y = '- log10 adjust p-value', title = 'tumor vs normal', color = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        panel.grid = element_blank(), 
        panel.background = element_rect(colour = 'black',fill = 'transparent')) +
  geom_vline(xintercept = c(-1,1), lty = 3, color = 'black') +
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-13, 7) + ylim(0,110)
p
ggsave("test_3\\TCGA_HNSC_volcano.pdf", p, width = 10, height = 8)

# 对显著表达的基因表达矩阵绘制热图，使用样本分组作为注释信息>>>
library(pheatmap)

expr_matrix = exp_gene[rownames(DEG_GENE),]
expr_matrix_scaled <- scale(expr_matrix, center = TRUE, scale = TRUE) #标准化
sample_Group = data.frame(group_list,row.names = colnames(expr_matrix_scaled)) #样本与分组对齐

head(expr_matrix_scaled)
max(expr_matrix_scaled)
min(expr_matrix_scaled)
length(expr_matrix_scaled)
length(expr_matrix_scaled[expr_matrix_scaled < -0.1])
expr_matrix_scaled_adj <- ifelse(expr_matrix_scaled >= 0.1, 0.1, expr_matrix_scaled)


p2 = pheatmap(
  expr_matrix_scaled_adj, 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = sample_Group,
  fontsize_row = 6,
  fontsize_col = 10,
  show_colnames = FALSE,
  show_rownames = FALSE,
)

library(ggplot2)
ggsave("test_3\\heatmap.png", p2, width = 10, height = 8,dpi = 300)
ggsave('test_3\\TCGA_HNSC_heatmap.pdf', p2, width = 10, height = 8)

