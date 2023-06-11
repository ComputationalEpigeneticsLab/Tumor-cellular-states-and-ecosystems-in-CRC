rm(list=ls())
library(clusterProfiler)
library(topGO)
# library(Rgraghviz)
library(pathview)
library(patchwork)
library(org.Hs.eg.db)
library(ggprism)
library(ggpubr)
library(ggplot2)
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}  
my_theme<-theme(panel.background=element_blank(),
                axis.text.y=element_text(color="black",size=10,family="Times",face="plain"),
                axis.text.x=element_text(colour="black",family="Times",size=10), 
                axis.title.y=element_text(family="Times",size = 10,face="plain"), 
                axis.title.x=element_text(family="Times",size = 10,face="plain"),
                plot.title = element_text(family="Times",size=12,face="bold",hjust = 0.5),
                panel.border = element_blank(),axis.line.x = element_line(colour = "black",size=0.5),
                axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                legend.text=element_text( family="Times", colour="black",size=10),
                legend.title=element_text( family="Times", colour="black",size=10))

road='F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\CRC的EcoTyper结果\\CRC_DiscoveryOutput\\'

addSmallLegend <- function(myPlot, pointSize = 2, textSize = 6, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
cells<-c('B.cells','CD4.T.cells','CD8.T.cells','Dendritic.cells','Endothelial.cells',
         'Epithelial.cells','Fibroblasts','Mast.cells','Monocytes.and.Macrophages','NK.cells',
         'PCs','PMNs')
setwd('F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos\\Hallmarker')
go_result_road<-'F:\\2study\\CRC_study\\data\\real_data\\发现细胞状态基因富集分析0.01\\GO_result\\'
for(i in 1:length(cells)){
  gene_info<-read.table(paste0(road,cells[i],'\\gene_info.txt'),header = T,as.is = T,sep = '\t')
  states<-unique(gene_info$State)
  for(j in 1:length(states)){
    state_gene<-gene_info[which(gene_info$State==states[j]),'Gene']
    state_gene<- bitr(state_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb='org.Hs.eg.db')
    state_gene<-na.omit(state_gene)
    state_go= enrichGO(gene = state_gene$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.01,
                       readable = T,
                       minGSSize = 15,maxGSSize = 500
    )
    go_gene<-as.data.frame(state_go@result)
    write.table(go_gene,paste0(go_result_road,cells[i],'_',states[j],'.txt'),quote = F,row.names = F,
                col.names = T,sep = '\t')
    go_gene$Count<-as.numeric(go_gene$Count)
    go_gene<-go_gene[order(go_gene$Count,decreasing =T),]
    go_gene<-head(go_gene,10)
    #library(patchwork)
    pdf(paste0(cells[i],'_',states[j],'.pdf'),width = 8,height = 3)
    p1=ggplot(data=go_gene,aes(x=Description,y=Count,fill=qvalue)) +
      geom_bar(stat="identity",alpha=0.8,position=position_dodge(0.75)) + coord_flip()+
      scale_fill_gradient(high="#4477AA",low="#BB4444")+
      # theme_prism()+
      scale_x_discrete(limits=rev(go_gene$Description)) +
      labs(x="",y="",title=paste0(cells[i],'_',states[j]))+
      scale_y_continuous(expand = c(0,0),breaks = integer_breaks())+
      my_theme
    p1=addSmallLegend(p1)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    print(p1+ plot_layout(ncol = 1))
    dev.off()
  }
  cat(i,sep = '\n')
}

