rm(list=ls())
library(ggplot2)
library(tidyverse)
#install.packages("see")
library(see)

dfc <- read.csv("恢复占比c.txt",stringsAsFactors = F,header = T,sep = "\t")
head(dfc)
pivot_longer(dfc,
             !Cell_Type,
             names_to = "Y",
             values_to = "Value") -> dfc.1
head(dfc.1)

library(ggplot2)
library(dplyr)
pdf("heat3.pdf",width = 7,height = 8)
dfc.1$Y<-factor(dfc.1$Y,
                levels = rev(c('B.cells','CD4.T.cells','CD8.T.cells','Dendritic.cells','Endothelial.cells',
                               'Epithelial.cells','Fibroblasts','Mast.cells','Monocytes.and.Macrophages',
                               'NK.cells','PCs','PMNs','EcoTyper','ALL')))


ggplot(dfc.1, aes(Cell_Type, Y)) + 
  geom_tile(aes(fill = Value),size=1,colour = "black")+
  scale_fill_gradient2(low = "#F5B4BF",high = "#EA5384")+
  geom_text(aes(label=Value),col ="black",size = 4)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5,vjust= 0.5,
                                   colour="black",family="Times",size=12),
        axis.text.y = element_text(size = 12,family="Times",face="plain"),
        
  )
dev.off() 
