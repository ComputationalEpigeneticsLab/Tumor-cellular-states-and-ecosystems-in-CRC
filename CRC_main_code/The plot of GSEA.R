


################################################################################################
######################one vs others #####################################
rm(list = ls())
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(doBy)###
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

discovery_exp<-read.csv("discovery_end_exp.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)

filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\cell_group_info")


cells<-c('B.cells','CD4.T.cells','CD8.T.cells','Dendritic.cells','Endothelial.cells',
         'Epithelial.cells','Fibroblasts','Mast.cells','Monocytes.and.Macrophages','NK.cells',
         'PCs','PMNs')


###
diff_result_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result\\"


Sys.time()
for (r in 1:length(cells)) {
  
  
 
  state<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\cell_group_info\\",
                           filename_select_group[r]),as.is=T,header = T,sep = "\t")
  cell_state<-state[,c(1,2)]
  
 
  gene<-discovery_exp[,which(colnames(discovery_exp) %in% cell_state$ID)]
  
  state_name<-unique(cell_state$State)
 
  for (m in 1:length(state_name)) {
    
    cell_state[which(cell_state$State==state_name[m]),3]<-c("group1")
    cell_state[!(cell_state$State==state_name[m]),3]<-c("group2")
    colnames(cell_state)<-c("sample","state","group")
    group<-cell_state[,c(1,3)]
    ad_state<-paste0(gsub('_(.*)','',filename_select_group[r]),"_",state_name[m])
    
    result<-c()
    
    for (n in 1:nrow(gene)) {
      
      gene_n <- data.frame(t(gene[n,]))
      gene_id <- names(gene_n)[1]
      names(gene_n)[1] <- 'gene'
      
      gene_n$sample <- rownames(gene_n)
      gene_n <- merge(gene_n, group, by = 'sample', all.x = TRUE)
      
      gene_n$group <- factor(gene_n$group)
      p_value <- wilcox.test(gene~group, gene_n)$p.value
      # if (!is.na(p_value) & p_value < 0.05) {
      stat <- summaryBy(gene~group, gene_n, FUN = c(mean, median))
      ###结果
      # result <- rbind(result, c(gene_id, as.character(stat[1,1]), stat[1,2], stat[1,3], 
      #                           as.character(stat[2,1]), stat[2,2], stat[2,3], p_value,state_name[m],ad_state))
      # 
      result <- rbind(result, c(gene_id, stat[1,2],stat[2,2],p_value,state_name[m],ad_state))
      # }
      
    }
    
    result <- data.frame(result)
    names(result) <- c('gene_id', 'mean1', 'mean2', 'p_value','state','cell_state')
    write.table(result,paste0(diff_result_road,cells[r],'_',state_name[m],'.txt'),
                quote = F,row.names = F,col.names = T,sep = '\t')
    
  }
  
  
}

Sys.time()




###
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")


filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_addpvalue\\"


for (i in 1:length(filename_select_group)) {
  
  
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  state_diff$p_adjust <- p.adjust(state_diff$p_value, method = 'BH')
  
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
  
}



###-logp
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_addpvalue")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logp\\"

###
for (i in 1:length(filename_select_group)) {
  
  
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_addpvalue\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    p_value<-state_diff[m,]$p_value
    logp<- -log10(p_value)
    state_diff[m,8]<-logp
    colnames(state_diff)[8]<-c("logp")
  }
  
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
  
}





###logFC
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logp")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logFC\\"

###
for (i in 1:length(filename_select_group)) {
  
  
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logp\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    mean1<-state_diff[m,]$mean1
    mean2<-state_diff[m,]$mean2
    FC<-mean1/mean2
    logFC<-log2(FC)
    state_diff[m,9]<-logFC
    colnames(state_diff)[9]<-c("logFC")
    
  }
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}





###sort
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")


filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logFC")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_sort\\"


for (i in 1:length(filename_select_group)) {
  
  #
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_logFC\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    logp<-state_diff[m,]$logp
    logFC<-state_diff[m,]$logFC
    sort<-logp*sign(logFC)
    state_diff[m,10]<-sort
    colnames(state_diff)[10]<-c("Sort")
    
  }
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}

##################################################################################################



##################################################################################
###########################GSEA###################################
####计算67个通路富集分析的GSEA
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

rm(list = ls())

##
eachstate_GSEA_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result\\"


setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

kegmt<-read.gmt("all_pathway.gmt")


###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_sort")

for (i in 1:length(filename_select_group)) {
  
  gene_list<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_sort\\",
                               filename_select_group[i]),as.is=T,header = T,sep = "\t")
  gene_list<-gene_list[,c(1,10,6)]
  gene_list[,2]=as.numeric(gene_list[,2])
  gene_list=gene_list[order(gene_list[,2],decreasing = T),]
  geneList<-gene_list[,2]
  names(geneList)<-as.character(gene_list[,1])
  KEGG<-GSEA(geneList,TERM2GENE = kegmt,
             minGSSize = 1,maxGSSize = 1000,
             pvalueCutoff = 1)
  KEGG<-as.data.frame(KEGG)
  cell_state<-unique(gene_list$cell_state)
  KEGG$cell_state<-cell_state
  KEGG_result<-KEGG
  
  KEGG_result_simply<-KEGG_result[,-c(1,10)]
  rownames(KEGG_result_simply)<-NULL
  KEGG_result_simply$Description<-tolower(KEGG_result_simply$Description)
  KEGG_result_simply$Description<-gsub("hallmark_","",KEGG_result_simply$Description)
  KEGG_result_simply$Description<-gsub("_"," ",KEGG_result_simply$Description)
  write.table(KEGG_result_simply,paste0(eachstate_GSEA_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
  
}


#################################################################################################
###################GSEA plot
###
rm(list = ls())
eachstate_GSEA_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_plot\\"
# setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\GSEA_result")


###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result")

for (i in 1:length(filename_select_group)) {
  
  GSEA_list<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result\\",
                               filename_select_group[i]),as.is=T,header = T,sep = "\t")
  GSEA_list$NES<-sort(GSEA_list$NES,decreasing = F)
  
  GSEA_list$Rank<-seq(1,nrow(GSEA_list),1)
  GSEA_list<-GSEA_list[,c(1,4,11)]
  
  pdf(paste0(gsub(".txt","",filename_select_group[i]),'.pdf'),width=12,height = 6)
  # GSEA_plot_data1<-GSEA_plot_data[which(GSEA_plot_data$NES>0),]
  ggplot(GSEA_list,
         mapping=aes(x=Rank,y=NES))+
    geom_point(size=2,shape=20)+
    # scale_color_manual(values = color)+
    scale_x_continuous(breaks = seq(1, 67, by = 1))+
    scale_y_continuous(breaks = seq(-2, 2.5, by = 0.1))+
    theme_bw()+theme(panel.grid = element_blank())
  
  
  dev.off()
  
  
}



###
###CD4.T.cells_S02
rm(list = ls())

setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result")

GSEA_list<-read.csv("CD4.T.cells_S02.txt",stringsAsFactors = F,header = T,sep = "\t")

###
GSEA_list=GSEA_list[order(GSEA_list$NES,decreasing = F),]


GSEA_list$Rank<-seq(1,nrow(GSEA_list),1)
GSEA_list<-GSEA_list[,c(1,4,6,11)]
# GSEA_list$Description<-tolower(GSEA_list$Description)
# GSEA_list$Description<-gsub("hallmark_","",GSEA_list$Description)
# GSEA_list$Description<-gsub("_"," ",GSEA_list$Description)
library(ggplot2)
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_plot")
pdf('CD4.T.cells_S02.pdf',width=14,height = 7)
# GSEA_plot_data1<-GSEA_plot_data[which(GSEA_plot_data$NES>0),]
p<-ggplot(GSEA_list,
          mapping=aes(x=Rank,y=NES,col=p.adjust))+
  geom_point(size=2,shape=19)+
  # scale_color_manual(values = color)+
  scale_x_continuous(breaks = seq(1, 67, by = 1))+
  scale_y_continuous(breaks = seq(-3.0, 2.8, by = 0.2))+
  theme_bw()+theme(panel.grid = element_blank())+
  # scale_color_continuous(low="",h)+
  scale_color_gradient(low = "#c60971",high = "#dd92bb")+
  labs(title = "CD4.T.cells S02")
# 
# # mygeneall<-GSEA_list[(nrow(GSEA_list)-9):nrow(GSEA_list),]
# library(gg.gap)
# gg.gap(plot = p,ylim = c(-3.0, 2.8),
#        segments = list(c(-0.8,0.4)),
#        tick_width = c(0.2,0.2,0.1),
#        rel_heights=c(0.3,0,0.5)
# )

mygeneall<-GSEA_list[GSEA_list$p.adjust<=0.001,]
library(ggrepel)
p+ geom_text_repel(
  data = GSEA_list[which(GSEA_list$Description %in% mygeneall$Description),],#res2_GSE134056[res2_GSE134056$pvalue<0.05&abs(res2_GSE134056$log2FoldChange)>0.25,],
  aes(label = Description),
  size =4,
  color = "black",#mycolor[res2[res2$padj<0.05&abs(res2$log2FoldChange)>1,'gene_type']],#
  segment.color = "black", show.legend = FALSE )

dev.off()




###Fibroblasts_S03
rm(list = ls())

setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result")

GSEA_list<-read.csv("Fibroblasts_S03.txt",stringsAsFactors = F,header = T,sep = "\t")

###
GSEA_list=GSEA_list[order(GSEA_list$NES,decreasing = F),]


GSEA_list$Rank<-seq(1,nrow(GSEA_list),1)
GSEA_list<-GSEA_list[,c(1,4,6,11)]
# GSEA_list$Description<-tolower(GSEA_list$Description)
# GSEA_list$Description<-gsub("hallmark_","",GSEA_list$Description)
# GSEA_list$Description<-gsub("_"," ",GSEA_list$Description)
library(ggplot2)
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_plot")
pdf('Fibroblasts_S03.pdf',width=14,height = 8)
# GSEA_plot_data1<-GSEA_plot_data[which(GSEA_plot_data$NES>0),]
p<-ggplot(GSEA_list,
          mapping=aes(x=Rank,y=NES,col=p.adjust))+
  geom_point(size=2,shape=19)+
  # scale_color_manual(values = color)+
  scale_x_continuous(breaks = seq(1, 67, by = 1))+
  scale_y_continuous(breaks = seq(-3.2, 1.8, by = 0.2))+
  theme_bw()+theme(panel.grid = element_blank())+
  # scale_color_continuous(low="",h)+
  scale_color_gradient(low = "#5310b2",high = "#8887b2")+
  labs(title = "Fibroblasts S03")
# 
# # mygeneall<-GSEA_list[(nrow(GSEA_list)-9):nrow(GSEA_list),]
# library(gg.gap)
# gg.gap(plot = p,ylim = c(-3.0, 2.8),
#        segments = list(c(-0.8,0.4)),
#        tick_width = c(0.2,0.2,0.1),
#        rel_heights=c(0.3,0,0.5)
# )

mygeneall<-GSEA_list[GSEA_list$p.adjust<=0.001,]
library(ggrepel)
p+ geom_text_repel(
  data = GSEA_list[which(GSEA_list$Description %in% mygeneall$Description),],#res2_GSE134056[res2_GSE134056$pvalue<0.05&abs(res2_GSE134056$log2FoldChange)>0.25,],
  aes(label = Description),
  size =4,
  color = "black",#mycolor[res2[res2$padj<0.05&abs(res2$log2FoldChange)>1,'gene_type']],#
  segment.color = "black", show.legend = FALSE )

##
# BiocManager::install("gg.gap")


dev.off()




####GSEA: top 5
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

rm(list = ls())

##
eachstate_GSEA_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result\\"



setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

kegmt<-read.gmt("all_pathway.gmt")
kegmt$term<-tolower(kegmt$term)
kegmt$term<-gsub("hallmark_","",kegmt$term)
kegmt$term<-gsub("_"," ",kegmt$term)


setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\diff_result_sort")

gene_list<-read.csv("Monocytes.and.Macrophages_S04.txt",stringsAsFactors = F,header = T,sep = "\t")
gene_list<-gene_list[,c(1,10,6)]
gene_list[,2]=as.numeric(gene_list[,2])
gene_list=gene_list[order(gene_list[,2],decreasing = T),]


geneList<-gene_list[,2]
names(geneList)<-as.character(gene_list[,1])
KEGG<-GSEA(geneList,TERM2GENE = kegmt,
           minGSSize = 1,maxGSSize = 1000,
           pvalueCutoff = 1)
library(ggthemes)
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_plot")
pdf("Mon_S04_top2.pdf",width = 8,height = 6)
color<-c(
  # '#1C9A35','#486D9C'
  # ,'#EC7A22',
  '#D53B40','#5FA69F'
)



gseaplot2(KEGG,c(39,44),rel_heights = c(0.8, 0.2, 0.6),
          base_size = 16,color = color,title = "Mon/Mac S04")
dev.off()





rm(list = ls())
# eachstate_GSEA_road<-"F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_plot\\"
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\GSEA_result")

filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result")


all_result<-c()
for (i in 1:length(filename_select_group)) {
  
  GSEA_list<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result\\",
                               filename_select_group[i]),as.is=T,header = T,sep = "\t")
  all_result<-rbind(all_result,GSEA_list)
  write.csv(all_result,"all_cell_states_GSEA_result.csv",row.names=F,col.names=T,quote=F)
  
}


GSEA_list<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\GSEA\\67_GSEA_result\\",
                             filename_select_group[i]),as.is=T,header = T,sep = "\t")







