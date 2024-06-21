library(TCGAbiolinks)
TCGAbiolinks:::getGDCprojects()$project_id
cancer_type="TCGA-HNSC"
clinical=GDCquery_clinic(project=cancer_type,type="clinical")

dim(clinical)
clinical[1:4,1:4]
head(colnames(clinical))

load("test_3/TCGA_HNSC_EXP.Rdata")
# ---临床与基因初始数据载入完毕---

#将第一列设置为行名
exp = gene_counts
group_list=ifelse(as.numeric(substr(colnames(exp),14,15)) < 10,'tumor','normal')
table(group_list)
exp_tumor = exp[,group_list == 'tumor']

# ---基因表达数据处理完毕---

meta = clinical
#同样以第一列作为列名
library(tibble)
meta <- column_to_rownames(meta,var = "submitter_id")
colnames(meta)
#筛选我们感兴趣的临床信息
meta=meta[,colnames(meta) %in% c("vital_status",
                                 "days_to_last_follow_up",
                                 "days_to_death",
                                 "gender",
                                 "age_at_index",
                                 "ajcc_pathologic_stage")]
dim(meta)
meta = meta[!(is.na(meta$days_to_last_follow_up) & is.na(meta$days_to_death)),]
meta = meta[meta$vital_status %in% c('Alive','Dead'),]

# ---临床数据处理完毕---

head(colnames(exp_tumor))
head(rownames(meta))
#调整、筛选临床样本信息数量与顺序与表达矩阵相同
colnames(exp_tumor) = gsub("\\.", "-", colnames(exp_tumor))
select_meta=meta[match(substr(colnames(exp_tumor),1,12),rownames(meta)),]

#1、计算生存时间
select_meta$days_to_death[is.na(select_meta$days_to_death)] <- 0   #缺失值标记为0
select_meta$days_to_last_follow_up[is.na(select_meta$days_to_last_follow_up)] <- 0
select_meta$days=as.numeric(select_meta[,'days_to_death'])+as.numeric(select_meta[,'days_to_last_follow_up'])
#时间以月份记，保留两位小数
select_meta$time=round(select_meta$days/30,2)

#2、根据生死定义活着是0，死的是1
select_meta$event=ifelse(select_meta$vital_status=='Alive',0,1)
table(select_meta$event)

#3 年龄分组(部分样本缺失，考虑可能的影响应该不大)
select_meta$age_at_index[is.na(select_meta$age_at_index)] <- 0
select_meta$age_at_index=as.numeric(select_meta$age_at_index)
select_meta$age_group=ifelse(select_meta$age_at_index>median(select_meta$age_at_index),'older','younger')
table(select_meta$age_group)

#4 癌症阶段
select_meta$ajcc_pathologic_stage[is.na(select_meta$ajcc_pathologic_stage)] <- 'unkown'
table(select_meta$ajcc_pathologic_stage)

#6 性别 gender
table(select_meta$gender)

dim(exp_tumor)
head(colnames(exp_tumor))
dim(select_meta)
head(rownames(select_meta))

exp_tumor = t(exp_tumor)
sur_data = cbind(select_meta,exp_tumor)

save(sur_data,file="test_3/tosur.RData")

# ---开始进行生存分析---
rm(list=ls())
load("test_3/tosur.RData")
library(survival)
library(survminer)
library(ggplot2)

#ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
#           risk.table =TRUE,pval =TRUE,
#           conf.int =TRUE,xlab ="Time in months", 
#           ggtheme =theme_light(), 
#           ncensor.plot = TRUE)

#肿瘤分期分组比较
sfit_tumor = survfit(Surv(time, event)~ajcc_pathologic_stage, data = sur_data)
tumor_curve = ggsurvplot(sfit_tumor,
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)# 手动保存吧……
tumor_curve

# 比较所有基因的log-rank p值
gene_rank = sur_data
gene_groups <- function(data, gene_name){
  group = ifelse(data[,gene_name]>median(data[,gene_name]),'high expression','low expression')
  return(group)

}
length(colnames(sur_data))# 11:13048
for (name in colnames(gene_rank)[11:13048]) {
  print(name)
  gene_rank[,name] = ifelse(gene_rank[,name]>median(gene_rank[,name]),'high expression','low expression')
}

gene_surv = with(gene_rank,Surv(time, event))
log_rank_p = c()
for (i in c(11:13048)){
  data_surv_diff = survdiff(gene_surv~gene_rank[,i],data = gene_rank)
  p.val = 1 - pchisq(data_surv_diff$chisq, length(data_surv_diff$n) - 1)
  log_rank_p =c(log_rank_p, p.val)
}
rankp_data = data.frame(colnames(gene_rank)[11:13048],log_rank_p)
rankp_data = rankp_data[order(rankp_data$log_rank_p),]
rankp_data = rankp_data[rankp_data$log_rank_p < 0.05,]
write.table(rankp_data,'test_3/TCGA_HNSC_prognostic_gene.txt',sep = '\t')
