rm(list=ls())

library(tximport)
library(DESeq2)

dir <- "./RNAseq/test/"
samples <- read.table(file.path(dir, "samples.txt"), header = F, stringsAsFactors = F)
colnames(samples)[1] <- 'run'
files <- file.path(dir, "", samples$run, "quant.sf")
all(file.exists(files))
for(i in seq(1,length(files))){ names(files)[i] <- strsplit(files,split="\\/")[[i]][7]}

mappings <- read.csv('./RNAseq/files/gencode.v29.mappings.tsv', sep = '\t', header = FALSE)
colnames(mappings) <- c('TXNAME', 'GENEID', 'type')
tx2gene <- mappings[c('TXNAME', 'GENEID')]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = T)

colnames(txi$abundance) <- gsub("/salmon","",samples$run)
colnames(txi$counts) <- gsub("/salmon","",samples$run)
colnames(txi$length) <- gsub("/salmon","",samples$run)

sampleTable <- data.frame(condition = factor(rep('A', ncol(txi$counts))))
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~1)
vst <- vst(dds, blind=TRUE)
vst <- assay(vst)

save(list = c("vst", "txi", ), file = "./RNAseq/files/strap_ready.RData")