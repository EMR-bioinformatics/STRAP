### DEG analysis for all three drugs

library(dplyr)
library(DESeq2)  
library(edgeR)    
library(qvalue)
library(purrr)

## loading data
load('RNASeq.RData')
metadata <- readRDS('metadata.rds')
table(metadata$Randomised.Medication, metadata$ACR.response.status)
# Etanercept   Rituximab Tocilizumab
#  67          72          69
txi <- STRAP_RNAseq_data$txi     # output from tximport with raw gene counts (gene count, not normalised)
## read in principal component table
pc <- read.delim('./PC_.txt', header = TRUE, sep = "\t", dec = ".")
pc <- subset(pc, select = -c(X, PC2, PC3))
## aligning metadata and data
identical(metadata$anonymised_Sample_Name, pc$SampleID)
metadata <- left_join(metadata, pc, by =c("anonymised_Sample_Name" = "SampleID"))
txi$abundance <- txi$abundance[, metadata$anonymised_Sample_Name]
txi$counts    <- txi$counts[,    metadata$anonymised_Sample_Name]
txi$length    <- txi$length[,    metadata$anonymised_Sample_Name] 
identical(metadata$anonymised_Sample_Name, colnames(txi$abundance))
parameter <- "ACR.response.status"  
drug.list <- c('Rituximab', "Tocilizumab", 'Etanercept')
summary_table <- lapply(drug.list, function(drug) {  
    # tidy & align data 
    metadata <- metadata %>% filter(
      Randomised.Medication == drug
    )
    metadata <- subset(metadata, !is.na(metadata[, parameter]))
    metadata[,parameter]  <-  as.factor(metadata[, parameter])  # make sure the column used for comparison is a factor 
    metadata[, parameter] <- droplevels(metadata[, parameter])
    txi$abundance <- txi$abundance[, metadata$anonymised_Sample_Name]
    txi$counts    <- txi$counts[,    metadata$anonymised_Sample_Name]
    txi$length    <- txi$length[,    metadata$anonymised_Sample_Name] 
    if(nlevels(metadata[, parameter]) == 2) {
      ########## DESeq2 pipeline ########## 
      dds <- DESeqDataSetFromTximport(txi, metadata, design= formula(paste0('~', pc, '+', parameter)))
      dds <- estimateSizeFactors(dds)       # find genes with low number of reads 
      idx <- rowSums( counts(dds, normalized=TRUE) >= 9 ) >= 18 ## filter out genes where there are less than 18 samples with normalized counts greater than or equal to 9.
      dds <- dds[idx,]
      dds <- DESeq(dds)
      res.deseq2 <- results(dds, contrast = c("ACR.response.status", "1", "0")) # 1 equals response
      ## adding q values (Storey and Tibshirani 2003)
      res.deseq2$qvalue <- NA
      res.deseq2$qvalue[!is.na(res.deseq2[, "padj"])] <- qvalue(res.deseq2[!is.na(res.deseq2[, "padj"]), "pvalue"])$qvalues
      res.deseq2 <- as.data.frame(res.deseq2)
      print(parameter)
      print(drug)
      print(nrow(subset(res.deseq2, res.deseq2$qvalue < 0.05)))
      #saveRDS(res.deseq2, file = paste0("deseq2.DEGs_", drug, "_filter918_PCadjust_", parameter, ".rds"))
    }
    return(dds)
  })
  
  return(summary_table_ACR20)
saveRDS(summary_table_ACR20, file="./ddsobjects_ACR20_filter918_PCadjust.rds")
