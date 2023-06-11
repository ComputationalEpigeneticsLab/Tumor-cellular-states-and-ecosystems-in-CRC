
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图")
data<-read.table("P_value.txt",stringsAsFactors = F,header = T)
data$Ecotypes<-factor(data$Ecotypes,levels = c("E7","E2","E5","E1","E3","E4","E6"))

pdf("survivial.pdf",width = 5,height = 3)
ggplot(data=data,
       mapping = aes(x=Ecotypes,y=-log10(number),fill=group),group=factor(1))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 12,face = "bold"),####
        axis.title.y = element_text(size = 12,face = "bold")
  )
dev.off()


###the abundance of cell types in seven ecotypes
rm(list = ls())
states_ecotype<-read.table("states_ecotype.txt",stringsAsFactors = F,header = T)
pivot_longer(states_ecotype,
             !cell_type,
             names_to = "Y",
             values_to = "Value") -> dfa
dfa$Y<-factor(dfa$Y,levels = c("E7","E2","E5","E1","E3","E4","E6"))

pdf("states_ecotype_abundance2.pdf",width = 5,height = 5)
ggplot(dfa, aes(Y,cell_type)) + 
  geom_tile(aes(fill = Value),size=1,colour = "black")+
  scale_fill_gradient2(low = "#FFFFFF",high = "#026C51")+
  # geom_text(aes(label=Value),col ="black",size = 4)+
  theme_minimal()+# 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),#
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5,vjust= 0.5,
                                   colour="black",family="Times",size=12),# 
        axis.text.y = element_text(size = 12,family="Times",face="plain"),
        
  )
dev.off() 


###
#the functional pathways of seven ecotypes
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图")
ecotype_gene<-read.table("ecotype_gene.txt",stringsAsFactors = F,header = T)
ecotype_gene<-ecotype_gene[,c(4,1)]

immuneGene<- readLines("immune.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))

mat<-matrix(ncol=5)
Ecotype<-unique(ecotype_gene[,1])

for (i in 1:length(Ecotype)) {
  
  for (j in 1:17) {
    
    n <- ecotype_gene[which(ecotype_gene[,1] == Ecotype[i]),2] #
    n_num <- length(n)
    m <- resGene[[j]] #
    m_num <- length(m)
    N_num <- 21876
    k <- intersect(n,m) #
    k_num <- length(k)
    pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
    #
    gene_pathway <- t(c(Ecotype[i], j, k_num, pvalues,paste(k,collapse=",")) )
    #
    write.table(gene_pathway, "gene_pathway.txt", row.names = F,
                col.names = F, sep = "\t", append=T, quote = F)
  }
  
}

gene_pathway2 <- read.table("gene_pathway.txt",header = F,sep = "\t",as.is = T,stringsAsFactors = F)

p<-p.adjust(gene_pathway2[,4],n=nrow(gene_pathway2),method = "fdr")###bonferroni
index1<-which(gene_pathway2[,4]<0.05)##
index2<-which(p<0.05)##
gene_pathway2_adjust<-gene_pathway2[index2,]
write.table(gene_pathway2_adjust,"gene_pathway2_adjust.txt",quote = F,sep = "\t",row.names = F)



##hallmarker
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图")
ecotype_gene<-read.table("ecotype_gene.txt",stringsAsFactors = F,header = T)
ecotype_gene<-ecotype_gene[,c(4,1)]
#
immuneGene<- readLines("Hallmarker.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))
#
mat<-matrix(ncol=5)
Ecotype<-unique(ecotype_gene[,1])

for (i in 1:length(Ecotype)) {
  
  for (j in 1:17) {
    
    n <- ecotype_gene[which(ecotype_gene[,1] == Ecotype[i]),2] #
    n_num <- length(n)
    m <- resGene[[j]] #
    m_num <- length(m)
    N_num <- 21876
    k <- intersect(n,m) #
    k_num <- length(k)
    pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
    #
    gene_pathway_hallmarker <- t(c(Ecotype[i], j, k_num, pvalues,paste(k,collapse=",")) )
    #
    write.table(gene_pathway_hallmarker, "gene_pathway_hallmarker.txt", row.names = F,
                col.names = F, sep = "\t", append=T, quote = F)
  }
  
}

gene_pathway_hallmarker2 <- read.table("gene_pathway_hallmarker.txt",header = F,sep = "\t",as.is = T,stringsAsFactors = F)

p<-p.adjust(gene_pathway_hallmarker2[,4],n=nrow(gene_pathway_hallmarker2),method = "fdr")###bonferroni矫正p值
index1<-which(gene_pathway_hallmarker2[,4]<0.05)#
index2<-which(p<0.05)##
gene_pathway_hallmarker2_adjust<-gene_pathway_hallmarker2[index2,]
write.table(gene_pathway_hallmarker2_adjust,"gene_pathway_hallmarker2_adjust.txt",quote = F,sep = "\t",row.names = F)


###GOBP
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图")
ecotype_gene<-read.table("ecotype_gene.txt",stringsAsFactors = F,header = T)
ecotype_gene<-ecotype_gene[,c(4,1)]

genelist<- bitr(ecotype_gene$Gene,fromType="SYMBOL",toType="ENTREZID",OrgDb='org.Hs.eg.db')
genelist2<-na.omit(genelist)
colnames(genelist2)<-c("Gene","ENTREZID")
library(clusterProfiler)
library(org.Hs.eg.db)
ecotype_GO_gene<-merge(ecotype_gene,genelist2,by="Gene")

ecotype<-c("E1","E2","E3","E4","E5","E6","E7")

#GOBP
go_result_road<-'F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图\\Ecotype_GOBP\\'
for(i in 1:length(ecotype)){
  ecotype_gene<-ecotype_GO_gene[which(ecotype_GO_gene$Ecotype==ecotype[i]),'ENTREZID']
  ecotype_go= enrichGO(gene = ecotype_gene,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable = T,
                       minGSSize = 15,maxGSSize = 500
  )
  go_gene<-as.data.frame(ecotype_go@result)
  write.table(go_gene,paste0(go_result_road,ecotype[i],'_GO_result','.txt'),quote = F,row.names = F,
              col.names = T,sep = '\t')
}

###
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\生存分析\\每个生态型生存分析p值柱形图")

pathway_info<-read.csv("enrichment_pathway2.txt",stringsAsFactors = F,header = T,sep = "\t")
# pathway_info<-pathway_info[,c(1,3,4,5)]
pathway_info$Ecotype<-factor(pathway_info$Ecotype,levels = c("E7","E2","E5","E1","E3","E4","E6"))

my_theme<-theme(axis.text.y=element_text(color="black",size=10,face="plain"),
                axis.text.x=element_text(colour="black",size=10), 
                # axis.title.y=element_text(size = 10,face="plain"), 
                # axis.title.x=element_text(size = 10,face="plain"),
                plot.title = element_text(size=12,face="bold",hjust = 0.5),
                axis.line.x = element_line(colour = "black",size=0.5),
                
                legend.text=element_text(colour="black",size=10),
                legend.title=element_text(colour="black",size=10))

pdf("ecotype_pathway2.pdf",width = 7,height =5)
ggplot(pathway_info,aes(Ecotype,pathway))+
  geom_point(aes(size=count,color=-log10(p.adjust.value)))+
  scale_colour_gradient(low = "#808080",high = "#000000")+
  theme_bw()+my_theme
dev.off()  





