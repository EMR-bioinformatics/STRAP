library(PCAtools)
library(Cairo)
plot.on.file = FALSE

# the below metadata was saved at line 101 of the "Figure4a_STRAP.clusters.R" script
metadata <- readRDS("../Data/STRAP/STRAP.metadata_clusters.RDS")
vst <- readRDS("../Data/STRAP/STRAP.vst.RDS")
vst <- vst[, rownames(metadata)]
vst <- vst[which(apply(vst, 1, var) != 0), ]

p <- PCAtools::pca(vst, metadata = metadata, scale = F, removeVar = 0.1)

# if(plot.on.file)(svg("Plots/STRAP.PCA.pathotypes.svg", width = 7.5, height = 6.5))
if(plot.on.file)(pdf("Plots/STRAP.PCA.pathotype.pdf", width = 7.5, height = 6.5))
biplot(p, x = "PC1", y = "PC2", 
       lab = NULL, axisLabSize = 22,
       colby = 'Pathotype',  # figure 4e-top
       colkey = c(Lymphoid = 'blue', Myeloid = 'red', Fibroid = "green"),
       # colby = 'cluster',  # figure 4e-bottom
       # colkey = c('cluster 1' = 'indianred2', 'cluster 2' = 'mediumseagreen',
       # 'cluster 3' = "steelblue"),
       ellipse = F,
       legendPosition = 'bottom',
       legendLabSize = 18,
       colLegendTitle = NULL
) + theme(axis.text = element_text(colour = "black")) 
if (plot.on.file) (dev.off())

