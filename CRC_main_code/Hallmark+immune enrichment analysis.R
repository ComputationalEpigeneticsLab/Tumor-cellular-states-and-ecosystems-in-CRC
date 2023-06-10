##hallmarker
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos")
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}  

#setwd("F:\\2study\\CRC_study\\data\\real_data\\发现状态标记基因富集")
road='F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\CRC的EcoTyper结果\\CRC_DiscoveryOutput\\'
#gene_info <- read.csv("B_gene_info.txt",stringsAsFactors = F,header = T,sep = "\t")
#rownames(gene_info) <- NULL
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


#Hallmarker
immuneGene<- readLines("Hallmarker.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))
ResMart<-c()
#setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos\\Hallmarker")
for (i in 1:length(cells)) {
  gene_info<-read.table(paste0(road,cells[i],'\\gene_info.txt'),header = T,as.is = T,sep = '\t')
  states<-unique(gene_info$State)
  ResMartMid3<-c()
  for (j in 1:length(states)) {
    states_gene<-gene_info[which(gene_info$State==states[j]),]$Gene
    ResMartMid2<-c()
    for (s in 1:50) {
      HallMarkerGenes<-resGene[[s]]
      n <- states_gene 
      n_num <- length(n)
      m <- HallMarkerGenes 
      m_num <- length(m)
      N_num <- 21876
      k <- intersect(n,m) 
      k_num <- length(k)
      pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
      #输出
      ResMartMid<-c(cells[i],states[j],names(resGene)[[s]],k_num,pvalues,paste(k,collapse=","))
      ResMart<-rbind(ResMart,ResMartMid)
      
    }
    
  }
  
}
colnames(ResMart)<-c('CellType','State','HallMarkerPathway','NumEnrich','pValue','Genes')
ResMart_lower<-as.data.frame(ResMart)
ResMart_lower$HallMarkerPathway<-tolower(ResMart_lower$HallMarkerPathway)



###
p<-p.adjust(ResMart_lower[,5],n=nrow(ResMart_lower),method = "fdr")###bonferroni
index1<-which(ResMart_lower[,5]<0.05)##
index2<-which(p<0.05)##
ResMart_adjust<-ResMart_lower[index2,]
write.table(ResMart_adjust,"ResMart_adjust_HallMarker_lower.txt",quote = F,sep = "\t",row.names = F)









##immune
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos")
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}  

#setwd("F:\\2study\\CRC_study\\data\\real_data\\发现状态标记基因富集")
road='F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\CRC的EcoTyper结果\\CRC_DiscoveryOutput\\'
#gene_info <- read.csv("B_gene_info.txt",stringsAsFactors = F,header = T,sep = "\t")
#rownames(gene_info) <- NULL
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


#immune
immuneGene<- readLines("Immune.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))
ResMart<-c()
#setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos\\Hallmarker")
for (i in 1:length(cells)) {
  gene_info<-read.table(paste0(road,cells[i],'\\gene_info.txt'),header = T,as.is = T,sep = '\t')
  states<-unique(gene_info$State)
  ResMartMid3<-c()
  for (j in 1:length(states)) {
    states_gene<-gene_info[which(gene_info$State==states[j]),]$Gene
    ResMartMid2<-c()
    for (s in 1:17) {
      HallMarkerGenes<-resGene[[s]]
      n <- states_gene 
      n_num <- length(n)
      m <- HallMarkerGenes 
      m_num <- length(m)
      N_num <- 21876
      k <- intersect(n,m) 
      k_num <- length(k)
      pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
      #输出
      ResMartMid<-c(cells[i],states[j],names(resGene)[[s]],k_num,pvalues,paste(k,collapse=","))
      ResMart<-rbind(ResMart,ResMartMid)
      
    }
    
  }
  
}
colnames(ResMart)<-c('CellType','State','immunePathway','NumEnrich','pValue','Genes')


###
p<-p.adjust(ResMart[,5],n=nrow(ResMart),method = "fdr")###bonferroni
index1<-which(ResMart[,5]<0.05)
index2<-which(p<0.05)
ResMart_adjust<-ResMart[index2,]
write.table(ResMart_adjust,"ResMart_adjust_immune_adjust.txt",quote = F,sep = "\t",row.names = F)







#####################################################################################################################
###The Circos of HallMarker+immune
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos\\富集分析")
##互作数据
states_pathway<-read.csv("HallMarker_immune_inter.txt",stringsAsFactors = F,header = T,sep = "\t")
# states_pathway<-read.csv("states_pathway.txt",stringsAsFactors = F,header = T,sep = "\t")

circos_inter<-states_pathway[,c(3,5,6)]
# circos_inter<-states_pathway
circos_inter<-circos_inter[which(circos_inter$HallMarkerPathway!=""),]



pathway_type<-read.csv("HallMarker_immune_type.txt",stringsAsFactors = F,header = T,sep = "\t")
pathway_type<-pathway_type[which(pathway_type$Node!=""),]


brand = c(structure(pathway_type$Type, names=pathway_type$Node))

brand = brand[!duplicated(names(brand))]
brand = brand[order(brand, names(brand))]
brand=factor(brand,levels = unique(brand)[c(1:7,10:14,8:9)])
brand = brand[order(brand)]

mycolor<-c("#E6AB02","#E7298A","#7570B3","#D95F02","#1B9E77","#CAB2D6",
           "#FDBF6F","#F7AEB3","#B2DF8A","#A6CEE3","#999999","#F780BF",
           "#D6372E","#FADD4B")

# 颜色
library(RColorBrewer)
library(randomcoloR)
# mycolor<-randomColor(15)#[c(1,2,4,5,6,7,8,9,11,12)]#[c(1:2,4:12)] [c(1,2,4,5,6,7,8,9,11,12)]
brand_color = structure(mycolor, names = levels(brand))
model_color = structure(rep(mycolor,as.numeric(table(brand))), names = names(brand))
pdf("Hallmarker_immune_circos2.pdf",width = 6,height = 6)
library(circlize)
circos.clear()
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(0.3, i-1), 4)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))
# 添加内部连线（有颜色）
chordDiagram(circos_inter, order = names(brand),
             grid.col = model_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02)
             ),
             annotationTrackHeight = mm_h(c(2, 2))
)
#
circos.track(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), sector.index, col = "black", cex =0.5, 
              facing = "reverse.clockwise", niceFacing = T)
}, bg.border = NA)
# 
for(b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(sector.index = model, track.index = 1, col = brand_color[b],
                   text = b, text.vjust = -1, niceFacing = TRUE
  )
}

circos.clear()
dev.off()



