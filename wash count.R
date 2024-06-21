library(TCGAbiolinks)
TCGAbiolinks:::getGDCprojects()$project_id
cancer_type="TCGA-HNSC"

query <- GDCquery(project = cancer_type,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query, directory = "bin", files.per.chunk = 50)
expdat <- GDCprepare(query = query,directory = "bin")
# >>>开始清洗>>>

library(AnnotationDbi)
library(SummarizedExperiment)
assays(expdat)
elementMetadata(expdat)$gene_name
gene_counts = assays(expdat)$fpkm_unstrand
rownames(gene_counts) = elementMetadata(expdat)$gene_name
length(table(rownames(gene_counts))) #59427

# 重复基因取平均
library(tibble)
gene_counts <- as_tibble(gene_counts, rownames = "Gene")
length(table(gene_counts$Gene)) #59427
library(dplyr)
merged_counts <- gene_counts %>%
  group_by(Gene) %>%
  summarise_all(mean)
length(table(merged_counts$Gene)) #59427
length(merged_counts$Gene)
# 过滤低表达量
gene_counts = merged_counts[!is.na(merged_counts$Gene),-1]
gene_counts = data.frame(gene_counts)
rownames(gene_counts) = merged_counts$Gene[!is.na(merged_counts$Gene)]

gene_counts = gene_counts[apply(gene_counts, 1, sum) >= length(colnames(gene_counts)),]
save(gene_counts, file = "test_3\\TCGA_HNSC_EXP.Rdata")
