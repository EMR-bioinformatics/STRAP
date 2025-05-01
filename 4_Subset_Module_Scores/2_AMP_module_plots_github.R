rm(list=ls())

useScDataGet <- "SCscores"

### Functions
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(circlize)
library(ComplexHeatmap)

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
  
  "SC-B1" = "IGHD+ CD27- naive",
  "SC-B2" = "IGHG3+ CD27+ memory",
  "SC-B3" = "Autoimmune associated",
  "SC-B4" = "Plasmablasts",
  
  "Cluster-0" = "TREM2low (0)",
  "Cluster-1" = "TREM2high (1)", 
  "Cluster-2" = "FOLR2+ ID2+ (2)",
  "Cluster-3" = "FOLR2high LYVE1+ (3)",
  "Cluster-8" = "FOLR2+ ICAM1+ (8)",
  "Cluster-4" = "HLAhigh CLEC10A+ (4)",
  "Cluster-7" = "HLAhigh ISG15+ (7)",
  "Cluster-5" =  "CD52+ S100A12+ (5)",
  "Cluster-6" = "CD52+ SSP1+ (6)")

### Prepare plots

load("SCscores.RData")
useScData <- get(useScDataGet)

all(colnames(useScData)==metadata$anonymised_Sample_Name)
all(metadata$anonymised_PatientID==Treatment$Patient.I.D.)
all(metadata$anonymised_PatientID==Labdata$Patient.I.D.)
all(metadata$anonymised_PatientID==myoutcome$Patient.I.D.)
myoutcome$anonymised_Sample_Name <- metadata$anonymised_Sample_Name
myoutcome$switched <- Treatment$switched

allResults <- list()
allResults_box <- list()

# Whether patient is responder or non-responder at that visit, according to ACR20 scores % change from baseline, using imputed (complete) data
r <-  "ACR.response.status"
table(myoutcome[[r]], useNA="always")
if(any(myoutcome[[r]] %in% c(1,0))){ myoutcome[[r]] <- ifelse(myoutcome[[r]]==1,"Responder","Non-responder") }

for(t in c("Tocilizumab", "Rituximab",  "Etanercept","Any")){
  
  if(t=="Any"){
    idxResponder <- which((myoutcome[[r]]=="Responder" |  myoutcome[[r]]==1) & Treatment$Randomised.Medication %in% c("Tocilizumab", "Rituximab",  "Etanercept"))
    idxNonResponder <- which((myoutcome[[r]]=="Non-responder" |  myoutcome[[r]]==0) & Treatment$Randomised.Medication %in% c("Tocilizumab", "Rituximab",  "Etanercept"))
  }else{
    idxResponder <- which((myoutcome[[r]]=="Responder" |  myoutcome[[r]]==1) & Treatment$Randomised.Medication==t)
    idxNonResponder <- which((myoutcome[[r]]=="Non-responder" |  myoutcome[[r]]==0) & Treatment$Randomised.Medication==t)
  }
  
  y <- useScData[,c(idxResponder, idxNonResponder)]  
  forModel <- c(myoutcome[[r]][idxResponder], myoutcome[[r]][idxNonResponder])
  forModel <- factor(forModel, levels = c("Non-responder", "Responder"))
    
  names(forModel) <- colnames(y)
  allResults[[t]] <- do.FC(y = y, forModel = forModel)
  y_melted <- reshape2::melt(y)
  allResults_box[[t]] <- cbind(y_melted, "Response"=forModel[y_melted$Var2]) 
  print(t)
  print(paste0(names(summary(forModel)), " ",summary(forModel)))
}

allResults_box_df <- do.call("rbind", allResults_box)
allResults_box_df$Treatment <- sapply(strsplit(rownames(allResults_box_df), split = "[.]"), function(x) x[1])
allResults_box_df$Var2 <- as.character(allResults_box_df$Var2) 
allResults_box_df$switched <- myoutcome$switched[match(allResults_box_df$Var2, myoutcome$anonymised_Sample_Name)]
myFC <- do.call("rbind", allResults)
rownames(myFC) <- gsub("[.]"," ", rownames(myFC))
rownames(myFC) <- gsub("_","-", rownames(myFC))

myFC$pval_symbl <- NA
myFC$pval_symbl[which(myFC$P.Value > 0.05)] <- ""
myFC$pval_symbl[myFC$P.Value < 0.05] <- "*"
myFC$pval_symbl[myFC$P.Value < 0.01] <- "**"
myFC$pval_symbl[myFC$P.Value < 0.001] <- "***"

myFC_ind_drug <- myFC[grep("Any|D2T",rownames(myFC), invert = T),]
myFC_any_ref <- myFC[grep("Any|D2T",rownames(myFC), invert = F),]

myFC_ind_drug$Treatment<- sapply(strsplit(rownames(myFC_ind_drug),split = " "), function(x) x[[1]])
myFC_ind_drug$subset <- myFC_ind_drug$celltype <- sapply(strsplit(rownames(myFC_ind_drug),split = " "), function(x) x[[2]])
myFC_ind_drug$subset <- factor(myFC_ind_drug$subset,
                               levels = c(grep("SC-F|SC-M",unique(myFC_ind_drug$subset), value = T),
                                          grep("SC-B|SC-T",unique(myFC_ind_drug$subset), value = T)))

myFC_ind_drug$celltype <- factor(myFC_ind_drug$celltype, 
                                 levels = c(grep("Fib|Mac",unique(myFC_ind_drug$celltype), value = T), 
                                            grep("B-cell|T-cell",unique(myFC_ind_drug$celltype), value = T)))

myFC_ind_drug$Treatment <- factor(myFC_ind_drug$Treatment, levels = c("Etanercept","Tocilizumab","Rituximab"))

colTreatment <- c("Etanercept"="#0d6db5", "Tocilizumab"="#c3280d", "Rituximab"="#2c9477")

p1 <- ggplot(data=myFC_ind_drug, aes(color=Treatment,
                                     x=reorder(subset, rev(as.numeric(as.factor(paste0(subset, celltype)))), decreasing = F),
                                     y=2^logFC, ymin=2^CI.L, ymax=2^CI.R, alpha=0.8)) + 
  geom_pointrange() +
  geom_text(size=8,
            label=myFC_ind_drug$pval_symbl, 
            nudge_x = -0.2, nudge_y = 0.1, 
            check_overlap = F
  ) + theme_bw(base_size = 18) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Fold change \n") +
  ggtitle("Responder vs. Non-responder") +
  facet_grid(. ~ Treatment, scales = "free", space = "free") + theme_classic2() + 
  scale_color_manual(values=colTreatment) + theme(axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"))
p1 <- p1 + theme(legend.position = "none")

p1
