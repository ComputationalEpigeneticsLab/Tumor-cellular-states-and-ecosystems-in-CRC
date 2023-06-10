

rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\状态和生态型的分析\\生态型中细胞状态的分布")

state_number<-read.csv("state_number2.txt",stringsAsFactors = F,header = T,sep = "\t")
library(tidyverse)



color<-c('B.cells'="#E6AB02", 'CD4.T.cells'="#E7298A", 'CD8.T.cells'="#7570B3", 
         'Dendritic.cells'="#D95F02",'Endothelial.cells'="#1B9E77",'Epithelial.cells'="#CAB2D6",
         'Fibroblasts'="#FDBF6F",'Mast.cells'="#F7AEB3",
         'Monocytes.and.Macrophages'="#B2DF8A",'NK.cells'="#A6CEE3",'PCs'="#999999",'PMNs'="#F780BF")
my_theme<-theme(panel.background=element_blank(),
                axis.text.y=element_text(color="black",size=10,family="Times",face="plain"),
                axis.text.x=element_text(colour="black",family="Times",size=10), 
                axis.title.y=element_text(family="Times",size = 10,face="plain"), 
                axis.title.x=element_text(family="Times",size = 10,face="plain"),
                plot.title = element_text(family="Times",size=12,face="bold",hjust = 0.5),
                panel.border = element_blank(),axis.line.x = element_line(colour = "black",size=0.5),
                axis.ticks.y = element_blank(),
                # axis.line.y = element_blank(),
                legend.text=element_text( family="Times", colour="black",size=10),
                legend.title=element_text( family="Times", colour="black",size=10))
pdf("state_percent.pdf",width = 7,height = 7)
state_number %>% 
  ggplot(aes(x = Ecotype, fill = Cell_type)) + 
  geom_bar(position = position_fill(),width = 0.8,alpha=0.9) + 
  scale_fill_manual(values = color)+ 
  theme_classic() + 
  labs(y = 'Percent') + 
  theme(legend.position="right")+scale_y_continuous(expand = c(0,0))+
  # theme(axis.text.x =element_text(size=12), 
  #       axis.text.y=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  my_theme
dev.off()

