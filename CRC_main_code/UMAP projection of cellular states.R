rm(list = ls())
library(umap)
setwd("F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\生态型欧氏距离")
data<-read.csv("CD4.txt",stringsAsFactors = F,header = T,sep = "\t")
group<-read.csv("CD4_state_assignment.txt",stringsAsFactors = F,header = T,sep = "\t")

iris.data<-data[,c(2:252)]
iris.labels<-group$State
iris.umap = umap(iris.data,metric = "euclidean")
data_plot<-data.frame(iris.umap$layout)
data_plot$label<-group$State
colnames(data_plot) <- c('X','Y','label') 

library(ggplot2)
ggplot(data_plot, aes(x=X, y=Y, colour=label)) + 
  geom_point(size=4)+
  scale_colour_manual(values = c("#E41A71","#379DB8"))+
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          panel.background = element_blank())








