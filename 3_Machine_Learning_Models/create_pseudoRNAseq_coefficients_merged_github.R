rm(list=ls())

# Here we use ETA,TOC and RTX merged counts to obtain global correlation coefficients

####################################
#### Convert nCounter to RNAseq ####
####################################

library(Metrics)
library(reshape2)
library(patchwork)

myDirPath <- "/3TR/Pipeline3TR/"
load(paste0(myDirPath, "/data/RNAseq_and_nCounter_merged_ETA_TOC_RTX_2023_12_04.RData"))

nanostring <- geneAnnot[geneAnnot$V1=="Endogenous",]
n_vst <- intersect(nanostring$V2, rownames(bulkvst))
bulkvst <- bulkvst[match(n_vst, rownames(bulkvst)),]
normDataBySyntheticOligo <- normDataBySyntheticOligo[match(n_vst, rownames(normDataBySyntheticOligo)),]

normDataBySyntheticOligo <- normDataBySyntheticOligo[,sampleAnnot$UniqueName2[sampleAnnot$UseInMachineLearning=="Yes"]]
bulkvst <- bulkvst[,sampleAnnot$UniqueName2[sampleAnnot$UseInMachineLearning=="Yes"]]

### Linear regression bulk vs nCounter

all(rownames(bulkvst)==rownames(normDataBySyntheticOligo))
all(colnames(bulkvst)==colnames(normDataBySyntheticOligo))

bulkvst_melted <- melt(bulkvst[,colnames(normDataBySyntheticOligo)])
normData_melted <- melt(normDataBySyntheticOligo)

nCounterEstimates <- list()
lm_res <- list()
myplots <- list()
lmdata_plotData <- list()

for(g in rownames(bulkvst)){
  
  if(g %in% bulkvst_melted$Var1){
    
    lmdata <- data.frame("bulk"=bulkvst_melted$value[bulkvst_melted$Var1==g],
                         "nCounter"=normData_melted$value[normData_melted$Var1==g])
    
    # predict bulk using nCounter
    lmFit <- lm(formula = bulk ~ nCounter, data = lmdata)
    lmFitSum <- summary(lmFit)
    lmFitSum$coefficients[,"Estimate"]
    y = (lmFitSum$coefficients[,"Estimate"][2] * lmdata$nCounter) + lmFitSum$coefficients[,"Estimate"][1]
    lmdata$RNAseq_from_nCounter <- y
    corcoef <- cor(y, lmdata$bulk, method = "pearson", use = "na.or.complete")
    myrmse = rmse(lmdata$bulk, y)
    
    pScatter_RNAseq_vs_nCounter1 <- ggplot(lmdata, aes(x=bulk, y=nCounter)) +  geom_point(shape=15,size=2) + 
      stat_smooth(method = "lm", col = "red") + ggtitle(paste0(g," | r: ", round(corcoef, 2))) + .px 
    
    pScatter_RNAseq_vs_nCounter2 <- ggplot(lmdata, aes(x=bulk, y=RNAseq_from_nCounter)) +  geom_point(shape=15,size=2) + 
      stat_smooth(method = "lm", col = "red") + ggtitle(paste0(g," | r: ", round(corcoef, 2))) + .px      
    
    pScatter_RNAseq_vs_nCounter3 <- ggplot(lmdata, aes(x=nCounter, y=RNAseq_from_nCounter)) +  geom_point(shape=15,size=2) + 
      stat_smooth(method = "lm", col = "red") + ggtitle(paste0(g," | r: ", round(corcoef, 2))) + .px      
    
    
    lmdata$corcoef <- corcoef
    lmdata_plotData[[g]] <- lmdata[,c("bulk","nCounter","corcoef")]
    myplots[[g]] <- wrap_plots(pScatter_RNAseq_vs_nCounter1, pScatter_RNAseq_vs_nCounter2, pScatter_RNAseq_vs_nCounter3)
    
    lm_res[[g]]$Beta_nCounter <- lmFitSum$coefficients[,"Estimate"][2]
    lm_res[[g]]$Intercept <- lmFitSum$coefficients[,"Estimate"][1]
    lm_res[[g]]$myrmse <- myrmse
    lm_res[[g]]$corcoef <- corcoef
    
    nCounterEstimates[[g]] <- y
  }else{
    print(g)
  }
}

nCounterEstimates_df <- as.data.frame(do.call("rbind",nCounterEstimates)) 
colnames(nCounterEstimates_df) <- colnames(normDataBySyntheticOligo)

lm_res_df <- as.data.frame(do.call("rbind",lm_res))
lm_res_df_Gene <- rownames(lm_res_df)
lm_res_df <- as.data.frame(apply(lm_res_df, 2, unlist))
rownames(lm_res_df) <- lm_res_df_Gene
lm_res_df_view <- as.data.frame(apply(lm_res_df, 2, function(x) {round(as.numeric(x), 2)}))
dimnames(lm_res_df_view) <- dimnames(lm_res_df)
lm_res_df_view[order(lm_res_df_view$corcoef, decreasing = T),]
all(bulkvst_melted$Var1==normData_melted$Var1)
all(bulkvst_melted$Var2==normData_melted$Var2)

lmdata <- data.frame("SampleID"=bulkvst_melted$Var2[bulkvst_melted$Var1 %in% rownames(normDataBySyntheticOligo)],
                     "Gene"=bulkvst_melted$Var1[bulkvst_melted$Var1 %in% rownames(normDataBySyntheticOligo)],
                     "bulk"=bulkvst_melted$value[bulkvst_melted$Var1 %in% rownames(normDataBySyntheticOligo)],
                     "nCounter"=normData_melted$value[normData_melted$Var1 %in% rownames(normDataBySyntheticOligo)])
lmdata$Run <- sampleAnnot$run[match(lmdata$SampleID, sampleAnnot$UniqueName2)]

lmFit <- lm(formula = bulk ~ nCounter, data = lmdata)
lmFitSum <- summary(lmFit)
lmFitSum$coefficients[,"Estimate"]
y = (lmFitSum$coefficients[,"Estimate"][2] * lmdata$nCounter) + lmFitSum$coefficients[,"Estimate"][1]
lmdata$y <- y
corcoef <- cor(y, lmdata$bulk, use = "na.or.complete", method = "pearson")
myrmse = Metrics::rmse(actual = lmdata$bulk[which(!is.na(y))], predicted = y[which(!is.na(y))])

# nCounterEstimates_df: pseudo RNAseq (generated using normDataBySyntheticOligo and lm_res_df)
# lm_res_df: Coefficients of linear regression (formula is bulkRNAseq = intercept + (beta * nCounter))

save(list=c("nCounterEstimates_df","lm_res_df"), 
     file=paste0(myDirPath,"/data/pseudoRNAseq_coefficients_merged_ETA_TOC_RTX_",format(Sys.Date(), format="%Y_%m_%d"),".RData"))

