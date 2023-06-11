rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\高表达的motif\\top5")
###
discovery_exp<-read.csv("discovery_end_exp.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
ecotype_sample<-read.csv("ecotype_assignment.txt",stringsAsFactors = F,header = T,sep = "\t")
ecotype_sample_exp<-discovery_exp[,which(colnames(discovery_exp) %in% ecotype_sample$ID)]

###TF
ecotype_TF<-read.csv("TF.txt",stringsAsFactors = F,header = T,sep = "\t")
TF<-unique(ecotype_TF$TF)
TF_expression<-ecotype_sample_exp[which(rownames(ecotype_sample_exp) %in% TF),]
write.table(TF_expression,"TF_expression_hot.txt",row.names = T,col.names = T,sep = "\t",quote=F)

###boxplot
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\高表达的motif\\top5")
##
TF_expression<-read.table("TF_expression.txt",stringsAsFactors = F,header = T)
ecotype_sample<-read.csv("ecotype_assignment.txt",stringsAsFactors = F,header = T,sep = "\t")
TF_expression_ecotype<-merge(ecotype_sample,TF_expression,by="ID")

one_ecotype_TF<-TF_expression_ecotype[,c(1,2,which(colnames(TF_expression_ecotype)=="POLR2A"))]

one_ecotype_TF[which(one_ecotype_TF$Ecotype=="CE6"),4]<-c("CE6")
one_ecotype_TF[!(one_ecotype_TF$Ecotype=="CE6"),4]<-c("Other CEs")
colnames(one_ecotype_TF)<-c("ID","Ecotype","TF","group")
# TF_expression_ecotype$V14
compare_means(TF~group,data = one_ecotype_TF,method = "t.test")

my_comparisons<-list(c("CE6","Other CEs"))

setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\高表达的motif\\top5\\TF_expression_boxplot")
pdf("CE6_POLR2A_point.pdf",width = 5,height = 5)
ggplot(one_ecotype_TF,
       aes(x=group,y=TF,fill=group))+ 
  stat_boxplot(geom = "errorbar",width=0.1)+ 
  # c('CE1'="#FF8B00",'Other CEs'="#0094C7")
  geom_boxplot(width=0.4,size=0.7,
               fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(color=group),width=.2,size=1)+
  
  # scale_fill_manual(values = c("red", "blue"))+  
  scale_color_manual(values=c("#D6372E","#0094C7"))+ 
  ggtitle("POLR2A_CE6")+ #
  theme_bw()+ #
  theme(legend.position="none", #
        axis.text.x=element_text(colour="black",family="Times",size=14), 
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #
        panel.grid.major = element_blank(), #
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  ylab("Expression")+xlab("Ecotype") #
dev.off()

###heatmap
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\高表达的motif\\top5")
#
TF_expression<-read.table("TF_expression_hot.txt",stringsAsFactors = F,header = T,row.names = 1)
#
ecotype_sample<-read.csv("ecotype_assignment.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
ecotype_sample$sample<-rownames(ecotype_sample)

ecotype_sample$Ecotype<-factor(ecotype_sample$Ecotype,levels = c('CE1','CE2','CE3','CE4','CE5','CE6','CE7'))
ecotype_sample<-ecotype_sample[order(ecotype_sample$Ecotype),]
ecotype_sample_fig<-data.frame('Ecotype'=ecotype_sample$Ecotype)
rownames(ecotype_sample_fig)<-rownames(ecotype_sample)

TF_expression<-TF_expression[,rownames(ecotype_sample_fig)]

gene_group<-read.csv("gene_group.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
gene_group$gene<-rownames(gene_group)

gene_group$gene<-factor(gene_group$gene,levels = c("IRF1","STAT1","IRF8","STAT2","IRF2","HSF1","HSF2",
                                                   "HSF4","MAFK","IRX5","FOXM1","SPIB","SPIC","SPI1","ELF5",
                                                   "GABPA","PHF8","TAF1","POLR2A","HCFC1","ELF1","SIN3A"))


gene_group<-gene_group[order(gene_group$gene),]
gene_group_fig<-data.frame('Ecotype'=gene_group$Ecotype)
rownames(gene_group_fig)<-rownames(gene_group)
TF_expression<-TF_expression[rownames(gene_group_fig),]

library(pheatmap)
pdf("TF_exp_heatmap.pdf",width = 8,height = 4.5)
bk=c(seq(-2,0,by=0.01),seq(0.01,2,by=0.01))
pheatmap(TF_expression, 
         # annotation_row=gene_group_fig, # 
         annotation_col=ecotype_sample_fig, # 
         show_colnames = F, #
         show_rownames=TRUE,  # 
         fontsize=8, # 
         # color = colorRampPalette(c('#4682B4','#FFFFFF','#D6372E'))(50), 
         color = c(colorRampPalette(colors = c("#4682B4","#FFFFFF"))(length(bk)/2),
                   colorRampPalette(colors = c("#FFFFFF","#D6372E"))(length(bk)/2)),
         annotation_legend=TRUE, #
         border_color=NA,  # 
         scale="row",  # 
         cluster_rows = F, #
         cluster_cols = FALSE, # 
         breaks=bk
         
)
dev.off()

