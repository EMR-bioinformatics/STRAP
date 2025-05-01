rm(list=ls())

library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MASS) 
library(sfsmisc)
library(ggridges)
library(ggstance)
options(Seurat.object.assay.version = "v3")

subsetLabels <- list(
  "SC-F1" = "CD34+ sublining",
  "SC-F2" = "HLA+ sublining",
  "SC-F3" = "DKK3+ sublining",
  "SC-F4" = "CD55+ lining",
  
  "SC-M1" = "IL1B+ pro-inflammatory",
  "SC-M2" = "NUPR1+",
  "SC-M3" = "C1QA+",
  "SC-M4" = "IFN-activated",
  
  "SC-T1" = "CCR7+ CD4+",
  "SC-T2" = "FOXP3+ Tregs",
  "SC-T3" = "PD-1+ Tph/Tfh",
  "SC-T4" = "GZMK+ CD8+",
  "SC-T5" = "GNLY+ GZMB+",
  "SC-T6" = "GZMK+/GZMB+",
  
  "SC-B1" = "IGHD+ CD270 naive",
  "SC-B2" = "IGHG3+ CD27- memory",
  "SC-B3" = "Autoimmune associated",
  "SC-B4" = "Plasmablasts",
  
  "Cluster-0" = "TREM2low (0)",
  "Cluster-1" = "TREM2high (1)", 
  "Cluster-2" = "FOLR2+   ID2+ (2)",
  "Cluster-3" = "FOLR2high LYVE1+ (3)",
  "Cluster-8" = "FOLR2+ ICAM1+ (8)",
  "Cluster-4" = "HLAhigh CLEC10A+ (4)",
  "Cluster-7" = "HLAhigh ISG15+ (7)",
  "Cluster-5" =  "CD52+  S100A12+ (5)",
  "Cluster-6" = "CD52+  SSP1+ (6)")

load("subsetModuleData_BulkRNAseq.RData")
rm(gexData, r4ra.meta, r4ra.switch)
load("/RNAseq/files/STRAP_db_v202202.RData")
gexData <- STRAP_db_v202202$RNAseq_data

######################
##  module scores  ##
#####################

strap_seurat <- CreateSeuratObject(counts = gexData, min.cells = 0, min.features = 0, project = "bulkRNAseq")
for(ct in c(unique(mygenes$cell.type))){
  mysubset <- mygenes[mygenes$cell.type==ct,]
  mysubset <- lapply(split(mysubset$gene, mysubset$clusterClean), function(x) x[1:5])  
  strap_seurat <- AddModuleScore(strap_seurat, features = mysubset, k=F, ctrl.size = 5, name = names(mysubset))  
}

names(strap_seurat@meta.data)[grep("SC_",names(strap_seurat@meta.data))] <- 
  sapply(grep("SC_",names(strap_seurat@meta.data), value = T), function(x) substr(x,1,nchar(x)-1))

names(strap_seurat@meta.data)[grep("Cluster",names(strap_seurat@meta.data))] <- 
  sapply(grep("Cluster",names(strap_seurat@meta.data), value = T), function(x) substr(x,1,nchar(x)-1))


# normalize the scores between 0-1
for(sset in unique(mygenes$clusterClean)){
  myscore <- as.vector(strap_seurat[[sset]][,1])
  mymin <- min(myscore)
  mymin <- ifelse(mymin > 0, -1 * mymin, abs(mymin))
  myscore <- myscore + mymin
  mymax <- max(myscore)
  myscore <- myscore/mymax
  strap_seurat[[sset]] <- myscore
}

SCscores <- t(strap_seurat@meta.data[,grep("SC_",names(strap_seurat@meta.data))])

save(list=c("SCscores"), file="/strap/result/SCscores.RData")
