###TCGA
library(survival)
library(survminer)
tfGene<-c('IRF1','STAT1','IRF8','STAT2','IRF2','HSF1','HSF2','HSF4','MAFK','IRX5'
          ,'FOXM1','SPIB','SPIC','SPI1','ELF5','GABPA','PHF8','TAF1','POLR2A'
          ,'HCFC1','ELF1','SIN3A')
setwd("F:\\2study\\CRC_study\\data\\real_data\\review2\\Lisi_CRC\\")
exprMat<-read.table('data\\TCGA_CRC_tumor2.txt',sep = '\t',as.is = T,
                    header = T,stringsAsFactors = F,check.names = F,row.names = 1)
exprMat<-log2(exprMat+1)
colnames(exprMat)<-gsub("\\.01","",colnames(exprMat))
surMat<-read.table('data\\TCGA_survival_data.txt',sep = '\t',as.is = T,
                   header = T,stringsAsFactors = F,check.names = F)
sampleUse<-intersect(colnames(exprMat),surMat$ID)

exprMatUse<-exprMat[tfGene,sampleUse]
gene_HR_res<-c()
for(i in 1:length(tfGene)){
  gene_exp<-as.data.frame(t(exprMatUse[tfGene[i],]))
  gene_exp$ID<-rownames(gene_exp)
  surMatGene<-merge(surMat,gene_exp,by='ID')
  colnames(surMatGene)<-c('ID','OS','OS.time','Expression')
  surMatGene$OS.time<-as.numeric(surMatGene$OS.time)
  res.cut <- surv_cutpoint(surMatGene, time = "OS.time", event = "OS",variables = "Expression")
  res.cat <- surv_categorize(res.cut)
  res.cat$Expression<-factor(res.cat$Expression,levels=c("low","high"))
  # fit <- survfit(Surv(OS.time,OS)~Expression,data = res.cat)
  # pdf(paste0('E:\\work\\project\\Lisi_CRC\\SuvivalPlotTCGA\\',tfGene[i],'_TCGASurPlot.pdf'),width=5,height=4,onefile = F)
  # p<-ggsurvplot(fit, data = res.cat,pval = TRUE,pval.size = 4,risk.table = F,pval.coord=c(100,10),
  #               xlab = "Overall Survival",  #x轴的label
  #               ylab = "Survival Proportion",font.subtitle = 12,font.tickslab = 12,
  #               font.x = 12,font.y = 12,fun = "pct",palette = c('#00BFC3','#F7766C'),
  #               submain=tfGene[i])
  # 
  # 
  # print(p)
  # dev.off()
  data.survdiff <- survdiff(Surv(OS.time,OS)~Expression,data = res.cat)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  ###合并
  gene_HR<-c(tfGene[i],p.val,HR,up95,low95)
  gene_HR_res<-rbind(gene_HR_res,gene_HR)
}

gene_HR_res<-as.data.frame(gene_HR_res)
colnames(gene_HR_res)<-c("TF","pvalue","HR","up95","low95")
gene_HR_res$fdr<-p.adjust(gene_HR_res$pvalue,n=length(gene_HR_res$pvalue),method = 'fdr')
#############################################################################################
gene_HR_res$pvalue<-as.numeric(gene_HR_res$pvalue)
gene_HR_res$HR<-as.numeric(gene_HR_res$HR)
gene_HR_res$low95<-as.numeric(gene_HR_res$low95)
gene_HR_res$up95<-as.numeric(gene_HR_res$up95)
#gene_HR_res<-gene_HR_res[order(gene_HR_res$HR),]
gene_HR_res$TF<-factor(gene_HR_res$TF,levels = rev(tfGene))

func_round<-function(v){return(round(v,3))}
gene_HR_res$pvalue<-func_round(gene_HR_res$pvalue)
write.csv(gene_HR_res,"TCGA.csv")
ggplot(data=gene_HR_res, aes(y=TF,color='red'),ylab=NA) + 
  geom_errorbar(data = gene_HR_res,aes(xmin = low95, xmax=up95, ylab=NA), #误差条表示95%的置信区间
                width=0.2,#误差条末端短横线的宽度
                position=position_dodge(0.9), 
                alpha = 1,xlim=c(0,3),
                size=1)  +theme(axis.text.y = element_blank())+
  theme_bw()+
  #绘制中位数为点图
  geom_point(data = gene_HR_res,aes(x=HR, y=TF),pch=19,position=position_dodge(0.9),size=4,ylab=NA,
             ylab=gene_HR_res$pvalue) +
  
  #scale_color_manual(values = c('#4393C3','#D6604D') )+ 
  #设置填充的颜色c("#56B4E9", "#E69F00") c("#1F78B4","#CB181D") c('#4393C3','#D6604D') 
  theme(
    panel.grid.major = element_blank(),   #不显示网格线
    panel.grid.minor = element_blank())+  #不显示网格线
  #xlab("Cancer")+ylab("log2(FPKM+1) (TRPA1)")+ #设置x轴和y轴的标题
  geom_vline(xintercept=1, linetype="dotted") 

