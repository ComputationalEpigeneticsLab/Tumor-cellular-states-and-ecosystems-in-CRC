rm(list = ls())
setwd("F:\\2study\\CRC_study\\data\\real_data\\IDWAS\\TCGA_drug_final")

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

GDSC2_Expr = readRDS('GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

 
testExpr=read.csv('TCGA_CRC_tumor2.txt',stringsAsFactors = F,header = T,sep = "\t")
rownames(testExpr)=testExpr$Gene
testExpr=testExpr %>% dplyr::select(-Gene)
testExpr=as.matrix(testExpr)
testExpr=log2(testExpr+1)
 
library(preprocessCore)
dir.create('GDSC2')
setwd('GDSC2')
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat   qn
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' ,
              rsq = F,
              cc =F,
              report_pc =F,
              pcr=F
)





