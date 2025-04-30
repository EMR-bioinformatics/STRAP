# Generate DEGs -----------------------------------------------------------
rm(list = ls())
library(DESeq2)
library(dplyr)
library(qvalue)
library(easylabel)
library(writexl)
library(tidyverse)

response.criteria <- "ACR.response.status.V7"
gene.set <- "protein_coding"
write_to_csv <- TRUE

metadata <- readRDS("STRAP.metadata.RDS")
txi <- readRDS("STRAP.txi.RDS") 
PCA <- readRDS("PCA_STRAP.RDS")

# tidy & align data 
metadata <- subset(metadata, metadata$Week == "0") # select baseline samples
metadata <- subset(metadata, !is.na(metadata[, response.criteria]))
metadata[, response.criteria]  <-  as.factor(metadata[, response.criteria])  
metadata[, response.criteria] <- droplevels(metadata[, response.criteria])
metadata <- metadata[which(rownames(metadata) %in% colnames(txi$counts)), ]
table(metadata[, response.criteria]) # 75 non-responder, 133 Responders 

if(gene.set == "protein_coding"){
  coding.genes <- read.csv("gene_type.txt", stringsAsFactors = FALSE)[, c(1,2)] %>%
    filter(Gene.type == "protein_coding") %>%
    pull(Gene.name)
  txi$abundance <- txi$abundance[which(rownames(txi$abundance) %in% coding.genes), rownames(metadata)]
  txi$counts    <- txi$counts[   which(rownames(txi$counts)    %in% coding.genes), rownames(metadata)]
  txi$length    <- txi$length[   which(rownames(txi$length)    %in% coding.genes), rownames(metadata)]
}

PCA <- PCA[rownames(metadata), ]
metadata$PC1 <- PCA$PC1

dds <- DESeqDataSetFromTximport(txi, metadata, design= formula(paste0('~PC1 +', response.criteria)))
dds <- DESeq(dds, parallel = F)
res.deseq2 <- results(dds)

# adding q values (Storey and Tibshirani 2003)
DESeq2.qvalues <- try(qvalue(p = res.deseq2[which(!is.na(res.deseq2$padj)), "pvalue"])[["qvalues"]])
if (class(DESeq2.qvalues) == "try-error") {
  DESeq2.qvalues <- qvalue(p = res.deseq2[which(!is.na(res.deseq2$padj)), "pvalue"],
                           pi0 = 1)[["qvalues"]]
}
res.deseq2$qvalue <- NA
res.deseq2[which(!is.na(res.deseq2$padj)), "qvalue"] <- DESeq2.qvalues
if(write_to_csv) write.csv(res.deseq2, file = paste0("DEGs.", response.criteria, ".PC1cov.csv"))

easyVolcano(res.deseq2, useQ = T, fdrcutoff = 0.05,
            Ltitle = "Non-responder",
            Rtitle = "Responder",
            width = 550, height = 450)

# Generate suppl table 
# to avoid a too long table we filter on p < 0.05
res.deseq2 <- res.deseq2 %>%
  filter(pvalue < 0.05) %>%
  arrange(qvalue) %>%
  rownames_to_column('Gene') %>%
  select("Gene", "baseMean", "log2FoldChange", "pvalue", "qvalue")
