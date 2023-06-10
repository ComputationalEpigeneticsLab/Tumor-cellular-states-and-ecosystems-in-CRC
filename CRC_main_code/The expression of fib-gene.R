rm(list=ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\Fib相关基因热图")

Fib_state_gene<-read.table("Fib_state_gene.txt",stringsAsFactors = F,header = T)



library(ggpubr)
compare_means(ITGBL1~State,data = Fib_state_gene,method = "t.test")
my_comparisons<-list(c("S03","S02"),c("S03","S01"),c("S03","S04"))

pdf("ITGBL1.pdf",width = 4,height = 4)
ggplot(Fib_state_gene,
       aes(x=State,y=ITGBL1,fill=State))+ 
  stat_boxplot(geom = "errorbar",width=0.1)+ 
  geom_boxplot(width=0.2,size=0.2,fill=c('S01'='#D4256B','S02'='#3B95AE','S03'='#5CA84B',"S04"='#754B9A'),outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  
  scale_color_manual(values=c("black","black","black"))+ 
  ggtitle("ITGBL1")+ 
  theme_bw()+ 
  theme(legend.position="none", 
        axis.text.x=element_text(colour="black",family="Times",size=14), 
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 14,face="plain"), 
        axis.title.x=element_text(family="Times",size = 14,face="plain"), 
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  ylab("Expression")+xlab("States") 

dev.off()



pdf("BGN.pdf",width = 4,height = 4)
ggplot(Fib_state_gene,
       aes(x=State,y=BGN,fill=State))+ 
  stat_boxplot(geom = "errorbar",width=0.1)+ 
  geom_boxplot(width=0.2,size=0.2,fill=c('S01'='#D4256B','S02'='#3B95AE','S03'='#5CA84B',"S04"='#754B9A'),outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  
  scale_color_manual(values=c("black","black","black"))+ 
  ggtitle("BGN")+ 
  theme_bw()+ 
  theme(legend.position="none", 
        axis.text.x=element_text(colour="black",family="Times",size=14), 
        axis.text.y=element_text(family="Times",size=14,face="plain"), 
        axis.title.y=element_text(family="Times",size = 14,face="plain"), 
        axis.title.x=element_text(family="Times",size = 14,face="plain"), 
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  ylab("Expression")+xlab("States") 

dev.off()


rm(list = ls())
setwd("/boot3/lisi/study2/CRC/scRNA-seq/Integrated")
object_integrated_res<-readRDS('object_integrated_res.rds')

pdf('BGN_sc.pdf',width=5,height=5)
VlnPlot(object_integrated_res,features='BGN',pt.size=0.0001,cols=colors,group.by="Cell_type")
dev.off()

pdf('ITGBL1_sc.pdf',width=5,height=5)
VlnPlot(object_integrated_res,features='ITGBL1',pt.size=0.0001,cols=colors,group.by="Cell_type")
dev.off()






