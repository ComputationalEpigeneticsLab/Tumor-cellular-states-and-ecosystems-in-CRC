##############################################
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\CRC数据发现和验证EcoTyper\\细胞类型中样本比例")

all_cell_assignment<-read.csv("all_cell_type_example.txt",stringsAsFactors = F,header = T,sep = "\t")
all_cell_assignment$number <- 1
all_cell_assignment$cell_type<-factor(all_cell_assignment$cell_type,
                                      levels = rev(c('B.cells','CD4.T.cells','CD8.T.cells','GSE29621','Dendritic.cells',
                                                     'Endothelial.cells','Epithelial.cells','Fibroblasts','Mast.cells',
                                                     'Mon/Mac','NK.cells','PCs','PMNs')))


library(MetBrewer)
color<-c("#E41A71","#379DB8","#5BAF4A","#7B4EA3","#FF7600","#FFC800","#A65328","#F780EC","#999999")

pdf("states_persent.pdf",width = 7,height = 8)
all_cell_assignment %>% 
  ggplot(aes(x = cell_type, fill = State)) + 
  geom_bar(position = position_fill(),width = 0.8) + 
  scale_fill_manual(values = color)+ 
  theme_classic() + 
  labs(y = 'Percent') + 
  theme(legend.position="right")+scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x =element_text(size=12), 
        axis.text.y=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  coord_flip()
dev.off()

