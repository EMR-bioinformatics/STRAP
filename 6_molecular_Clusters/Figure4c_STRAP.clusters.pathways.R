rm(list = ls())
library(enrichR)
library(dplyr)
library(writexl)
library(Cairo)
source("../plot_pathways.R")
plot.on.file <- FALSE
cluster <- "1" # 1, 2, or 3

cluster.genes <- read.csv(paste0("STRAP.genes.cluster", cluster, ".csv"))[,1]

pathways <- enrichr(cluster.genes, "Reactome_2022")
pathways[["Reactome_2022"]]$Term <- unlist(lapply(strsplit(pathways[["Reactome_2022"]]$Term, "R-HSA"), `[[`, 1))
pathways <- pathways$Reactome_2022 %>%
  filter(Adjusted.P.value < 0.01)
write_xlsx(pathways, paste0("STRAP.pathways.cluster", cluster, ".xlsx"))

# Merged pathways plot ------------------------------------------------------
# In the previous section I generated the pathways associated to 
# each cluster using enrichR. I then exported an excel file for each of them.
# Together with Myles, the most interesting pathways were selected for plotting
# (see the "Interesting" columns)

cluster1.pathways <- read.csv("STRAP.pathways.cluster1.csv") %>% 
  filter(Interesting == "y") %>%
  mutate(Cluster = "Cluster 1")
cluster2.pathways <- read.csv("STRAP.pathways.cluster2.csv") %>% 
  filter(Interesting == "y") %>%
  mutate(Cluster = "Cluster 2") 
cluster3.pathways <- read.csv("STRAP.pathways.cluster3.csv") %>% 
  filter(Interesting == "y") %>%
  mutate(Cluster = "Cluster 3")  

# make sure we remove any leading or trailing white space
cluster1.pathways$Term <- trimws(cluster1.pathways$Term)
cluster2.pathways$Term <- trimws(cluster2.pathways$Term)
cluster3.pathways$Term <- trimws(cluster3.pathways$Term)

# shorten one term to fit the space
cluster1.pathways$Term[which(
  cluster1.pathways$Term ==  "L13a-mediated Translational Silencing Of Ceruloplasmin Expression")] <-
  "L13a-mediated Translational Silencing Of Ceruloplasmin"
cluster3.pathways$Term <- gsub("Fcgamma", "Fc-\u03B3", cluster3.pathways$Term)

all.pathways <- rbind(cluster1.pathways, cluster2.pathways, cluster3.pathways)

# invert row order for the plot 
all.pathways <- all.pathways[nrow(all.pathways):1, ]

CairoPDF("plots/STRAP.pathways.merged.pdf", width = 6, height = 8) # Cairo is needed to plot the gamma symbol (γ). 
plot.terms.dots.categorical(all.pathways, FDR.correction = T,
                            top.n = 30, # max number of pathways per category (cluster)
                            term.size = 13,
                            order.by.significance = FALSE,
                            color.criteria = 'Cluster',
                            line.color = 'darkgray',
                            dot.size = 15) + 
  scale_color_manual(values = list('Cluster 1' = "indianred2",
                                   'Cluster 2'  = "mediumseagreen",
                                   'Cluster 3'   = "steelblue"))  +
  theme(
    legend.position = "bottom",
    legend.justification = c(1.5, 0),
    legend.box.spacing = unit(0.04, "cm")
  )
if(plot.on.file) {dev.off()}