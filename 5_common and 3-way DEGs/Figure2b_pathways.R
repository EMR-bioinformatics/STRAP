rm(list = ls())
library(enrichR)
library(dplyr)
source("../plot_pathways.R")
plot.on.file <- TRUE

responder.DEGs <- read.csv(paste0("DEGs.ACR.response.status.V7.PC1cov.csv"), row.names = 1) %>% 
  filter(qvalue < 0.05) %>%
  filter(log2FoldChange > 0) %>%
  rownames()

non.responder.DEGs <- read.csv(paste0("DEGs.ACR.response.status.V7.PC1cov.csv"), row.names = 1) %>% 
  filter(qvalue < 0.05) %>%
  filter(log2FoldChange < 0) %>%
  rownames()

# enrichR
responder.pathways <- enrichr(responder.DEGs, "Reactome_2022")[[1]]
non.responder.pathways <- enrichr(non.responder.DEGs, "Reactome_2022")[[1]]

# Responders/Non-Responders Pathways 
responder.pathways2 <-  responder.pathways %>% 
  filter(P.value < 0.05) %>% 
  arrange(P.value) %>%
  select(-c(Old.P.value, Old.Adjusted.P.value))

non.responder.pathways2 <-  non.responder.pathways %>% 
  filter(P.value < 0.05) %>% 
  arrange(P.value) %>%
  select(-c(Old.P.value, Old.Adjusted.P.value))

responder.pathways2$P.value          <- round(responder.pathways2$P.value, digits = 3)
responder.pathways2$Adjusted.P.value <- round(responder.pathways2$Adjusted.P.value, digits = 3)
responder.pathways2$Odds.Ratio       <- round(responder.pathways2$Odds.Ratio, digits = 2)
responder.pathways2$Combined.Score   <- round(responder.pathways2$Combined.Score, digits = 2)
non.responder.pathways2$P.value          <- round(non.responder.pathways2$P.value, digits = 3)
non.responder.pathways2$Adjusted.P.value <- round(non.responder.pathways2$Adjusted.P.value, digits = 3)
non.responder.pathways2$Odds.Ratio       <- round(non.responder.pathways2$Odds.Ratio, digits = 2)
non.responder.pathways2$Combined.Score   <- round(non.responder.pathways2$Combined.Score, digits = 2)

# Supplementary Data 3
datalist <- list("DEGs" = DEGs,
  "Responders_Pathways" = responder.pathways2,
  "NonResponders_Pathways"  = non.responder.pathways2)
names(datalist) <- substr(names(datalist), start = 1, stop = 29)
write_xlsx(datalist, "Supplementary Data 3.xlsx")

responder.pathways$Term <- unlist(lapply(strsplit(responder.pathways$Term, "R-HSA"), `[[`, 1))
non.responder.pathways$Term <- unlist(lapply(strsplit(non.responder.pathways$Term, "R-HSA"), `[[`, 1))

# Abbreviating some long pathway name
responder.pathways$Term[2] <- "Antigen activates B Cell Receptor generating second messengers "
responder.pathways$Term[8] <- "Interactions between a Lymphoid and a non-Lymphoid cell "

# let's keep the most interesting ones for the figure
non.responder.pathways <- non.responder.pathways[-c(1,2,7,9), ]

# if(plot.on.file) {svg(paste("merged.pathways", response.criteria,
#                             "svg", sep = "."), width = 7, height = 7)}
if(plot.on.file) {pdf(paste("_merged.pathways", response.criteria,
                            "pdf", sep = "."), width = 7, height = 7)}
plot.merged.dots(pathway.table.left = responder.pathways,
                 pathway.table.right  = non.responder.pathways,
                 left.label =  "Responders" ,
                 right.label = "Non Responders",
                 pathway.tool = "enrichR",
                 top.n.left = 8,
                 top.n.right = 6,
                 term.size = 13,
                 title.position = 1,
                 plot.title = "",
                 FDR.correction = TRUE) + 
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10, 15),
                     labels = c("10", "5", "0", "5", "10", "15")) +
  theme(legend.position = "bottom",
        legend.justification = c(1.5, 0),
        legend.box.spacing = unit(0.04, "cm"),
        legend.text = element_text(colour = 'black', family = "sans", size = 10))
if(plot.on.file) {dev.off()}

