
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
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因")
###
discovery_exp<-read.csv("discovery_end_exp.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
ecotype_sample<-read.csv("ecotype_assignment.txt",stringsAsFactors = F,header = T,sep = "\t")

ecotype_name<-unique(ecotype_sample$Ecotype)

diff_result_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_gene\\"
for (m in 1:length(ecotype_name)) {
  
  gene<-discovery_exp[,which(colnames(discovery_exp) %in% ecotype_sample[,1])]

  ecotype_sample_group<-ecotype_sample
  ecotype_sample_group[which(ecotype_sample_group$Ecotype==ecotype_name[m]),3]<-c("group1")
  ecotype_sample_group[!(ecotype_sample_group$Ecotype==ecotype_name[m]),3]<-c("group2")
  colnames(ecotype_sample_group)<-c("sample","state","group")
  group<-ecotype_sample_group[,c(1,3)]
 
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
 
    result <- rbind(result, c(gene_id, stat[1,2],stat[2,2],p_value,ecotype_name[m]))
    # }
    
  }
  
  result <- data.frame(result)
  names(result) <- c('gene_id', 'mean1', 'mean2', 'p_value','Ecotype')
  write.table(result,paste0(diff_result_road,ecotype_name[m],'.txt'),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}


###
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_gene")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_addpvalue\\"

###
for (i in 1:length(filename_select_group)) {
  
  #
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_gene\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  state_diff$p_adjust <- p.adjust(state_diff$p_value, method = 'fdr')
  
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
  
}


###-logp
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_addpvalue")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logp\\"

###
for (i in 1:length(filename_select_group)) {
  
  #
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_addpvalue\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    p_value<-state_diff[m,]$p_value
    logp<- -log10(p_value)
    state_diff[m,7]<-logp
    colnames(state_diff)[7]<-c("logp")
  }
  
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
  
}


###logFC
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logp")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logFC\\"

for (i in 1:length(filename_select_group)) {
  
  #
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logp\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    mean1<-state_diff[m,]$mean1
    mean2<-state_diff[m,]$mean2
    FC<-mean1/mean2
    logFC<-log2(FC)
    state_diff[m,8]<-logFC
    colnames(state_diff)[8]<-c("logFC")
    
  }
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}


###sort
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logFC")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_sort\\"

for (i in 1:length(filename_select_group)) {
  
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_logFC\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  for (m in 1:nrow(state_diff)) {
    
    logp<-state_diff[m,]$logp
    logFC<-state_diff[m,]$logFC
    sort<-logp*sign(logFC)
    state_diff[m,9]<-sort
    colnames(state_diff)[9]<-c("Sort")
    
  }
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}
###
##############################################################################################


###
rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\GSEA")

###
filename_select_group <- list.files("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_sort")

diff_addpvalue_road<-"F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\ecotype_high_express_sig_FDR\\"


for (i in 1:length(filename_select_group)) {
  
  #
  state_diff<-read.table(paste0("F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\diff_result_sort\\",
                                filename_select_group[i]),as.is=T,header = T,sep = "\t")
  
  state_diff<-state_diff[which(state_diff$Sort>3),] ###p<0.001
  state_diff<-state_diff[which(state_diff$p_adjust<0.001),] ###FDR<0.05
  
  
  write.table(state_diff,paste0(diff_addpvalue_road,gsub(".txt","",filename_select_group[i]),".txt"),
              quote = F,row.names = F,col.names = T,sep = '\t')
  
}


rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif")

library(AUCell)
library(RcisTarget)


###############################################################################
###
rm(list = ls())
# setwd("/boot3/lisi/study2/CRC/motif")
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif")
###
txtFile <- file.path('F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\ecotype_high_express_sig_FDR\\E7.txt')
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

#motif-TF 
data(motifAnnotations_hgnc) # human TFs (for motif collection 9)

#
motifRankings <- importRankings("hg19-500bp-upstream-7species.mc9nr.feather")

# 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings,
                      aucMaxRank = 0.05 * ncol(motifRankings),nCores = 1)



# 2. Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                           motifAnnot=motifAnnotations_hgnc)


# 3. Identify significant genes for each motif
# (i.e. genes from the gene set in the top of the ranking)
# Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="aprox")

###motif
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\ecotype_motif_result")
pdf("E1_AUC_hist.pdf",width = 5,height = 5)
auc <- getAUC(motifs_AUC)[1,]
hist(auc, main="hypoxia", xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred",ylim = c(0,1400))
nes3 <- (3*sd(auc)) + mean(auc) 
abline(v=nes3, col="red")
dev.off()


###motif
data(motifAnnotations_hgnc)
motifAnnotations_hgnc
cg=auc[auc>nes3]
names(cg)
cgmotif=motifAnnotations_hgnc[match(names(cg),motifAnnotations_hgnc$motif),]
cgmotif=na.omit(cgmotif)

###
# motifEnrichmentTable_wGenes
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
library(DT)
datatable(motifEnrichmentTable_wGenes_wLogo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

###
anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

anotatedTfs$hypoxia
signifMotifNames <- motifEnrichmentTable$motif[1:3]

###
pdf("E1_AUC.pdf",width = 5,height = 5)
incidenceMatrix <- getSignificantGenes(geneLists$hypoxia, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix
dev.off()

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
# BiocManager::install("visNetwork")
library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))

pdf("E1_network.pdf",width = 15,height = 15)
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)
dev.off()


###
rm(list = ls())
# setwd("/boot3/lisi/study2/CRC/motif")
setwd("F:\\2study\\CRC_study\\data\\real_data\\motif")
###
txtFile <- file.path('F:\\2study\\CRC_study\\data\\real_data\\motif\\EcoType差异表达基因\\ecotype_high_express_sig_FDR\\E1.txt')
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

# txtFile <- file.path('/boot3/lisi/study2/CRC/motif/ecotype_high_express_sig/E1.txt')
# geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

motifRankings <- importRankings("hg19-500bp-upstream-7species.mc9nr.feather")

#
data(motifAnnotations_hgnc) # human TFs (for motif collection 9)
# motifAnnotations_hgnc[199:202,]
# motifAnnotation <- motifAnnotations_hgnc


motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc,
                                         nesThreshold=3, geneErnMethod="aprox", nCores=1,
                                         aucMaxRank = 0.05 * ncol(motifRankings),
                                         geneErnMmaxRank = 5000
)

# motifEnrichmentTable_wGenes2=na.omit(motifEnrichmentTable_wGenes)

setwd("F:\\2study\\CRC_study\\data\\real_data\\motif\\ecotype_motif_result")
write.table(motifEnrichmentTable_wGenes,"E7_motif.txt",quote = F,sep = "\t",row.names = F)

motifEnrichmentTable_wGenes_logo <- addLogo(motifEnrichmentTable_wGenes)
write.table(motifEnrichmentTable_wGenes_logo,"E7_motif_logo.txt",quote = F,sep = "\t",row.names = F)
library(DT)

###
datatable(motifEnrichmentTable_wGenes_logo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))
