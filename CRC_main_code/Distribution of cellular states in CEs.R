
rm(list=ls())
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(reshape)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("randomcoloR")

setwd("F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\状态和生态型的分析\\生态型中细胞状态的分布")

state_ecotype<-read.table("ecotypes.txt",stringsAsFactors=F,header=T)
state_ecotype<-state_ecotype[,c(1,5)]
state_ecotype$link<-1

state_ecotype<-reshape::melt(state_ecotype,id='link')
variable<-summary(state_ecotype$variable)
state_ecotype$flow<-rep(1:variable[1],length(variable))
head(state_ecotype)
write.table(virus_pathway,"virus_pathway.txt",quote = F,sep = "\t")

####sankey
library(ggalluvial)
library(RColorBrewer)
library(randomcoloR)

col_sankey<-c('B.cells'="#E6AB02", 'CD4.T.cells'="#E7298A", 'CD8.T.cells'="#7570B3", 
              'Dendritic.cells'="#D95F02",'Endothelial.cells'="#1B9E77",'Epithelial.cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast.cells'="#F7AEB3",
              'Monocytes.and.Macrophages'="#B2DF8A",'NK.cells'="#A6CEE3",'PCs'="#999999",'PMNs'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3",'E7'="#B7D3E5")

pdf("state_ecotype_sankey2.pdf",width = 6,height = 8)
ggplot(state_ecotype, aes(x = variable, y = link,
                          stratum = value, alluvium = flow, fill = value)) +
  geom_stratum() + 
  geom_flow(aes.flow = 'forward') + 
  scale_fill_manual(values = col_sankey[rank(1:44)]) + 
  guides(fill=FALSE)+
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) + 
  labs(x = '', y = '') + 
  theme(legend.position = 'none', panel.background = element_blank(),
        line = element_blank(), axis.text.y = element_blank())+
  scale_x_discrete(limits = c('CellType', 'Ecotype'))

dev.off()  











