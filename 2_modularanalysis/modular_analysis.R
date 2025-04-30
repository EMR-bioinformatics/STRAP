### Modular analysis for ACR20 using object with PC1 adjusted DEG results
library(DESeq2)
library(edgeR)
library(dplyr)
library(qusage)
library(qvalue)
source('load_mods_qmod.R') # load modules
source('qmodS.R') # load qmod function
#load data
summary_table_ACR20 <- readRDS("ddsobjects_ACR20_filter918_PCadjust.rds")
drug <- 'Etanercept' # select each drug individually, Etanercept, Rituximab, Tocilizumab
data <- summary_table_ACR20$Etanercept


mod.list <- c('reduced_wgcna','LI_reduced')
# all three drugs run this way
qusageres_list_ACR20_ETA<- lapply(mod.list, function(module){
  c2.indices <- loadmodsforqusage(load = module)
  qusageres <- qmod(fit=data, filter= c(1,1), coef=3, geneSets = c2.indices) #coef = 3 for analysis with covariate, PC1 has been added as covariate
  print(drug)
  print(module)
  print(nrow(subset(qusageres, qusageres$qval < 0.05)))
  print(nrow(subset(qusageres, qusageres$p.Value < 0.05)))
  
  return(qusageres)
})

names(qusageres_list_ACR20_ETA) <- unlist(mod.list)
saveRDS(qusageres_list_ACR20_ETA, file = 'ACR20_918_ETA_PCadjust_qusageres.rds')


