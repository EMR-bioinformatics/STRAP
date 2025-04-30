rm(list = ls())
library(dplyr)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# load & filter data ---------------------------------------------------------------
metadata <- readRDS("STRAP.metadata.RDS") %>% 
  filter(Week == 0) # baseline only
counts <- readRDS("STRAP.txi.RDS")$counts

coding.genes <- read.csv("gene_type.txt", stringsAsFactors = FALSE)[, c(1,2)] %>%
  filter(Gene.type == "protein_coding") %>%
  pull(Gene.name)
coding.genes <- coding.genes[coding.genes %in% rownames(counts)]

counts <- counts[which(rownames(counts) %in% coding.genes), rownames(metadata)]
keep <- filterByExpr(counts, min.count = 20, min.total.count = 500000)
table(keep)
counts <- counts[keep, ] # 14435 genes, 3411 for more stringent filterbyExpr  
nrow(counts)
mode(counts) <- "integer"

vst <- vst(counts, blind=TRUE) # do the vst transform, unsupervised

# Heatmap -----------------------------------------------------------------
plot.on.file = FALSE

heatmap.df <- t(scale(t(vst)))
heatmap.df <- heatmap.df[which(complete.cases(heatmap.df)),]

set.seed(101)

ha = HeatmapAnnotation(
      "CD3" = metadata$CD3,
      "CD20"  = metadata$CD20 ,
      "CD68L" = metadata$CD68L,
      "CD68SL" = metadata$CD68SL,
      "CD138"  = metadata$CD138,
      "Cell type" = metadata$Cell.type,
      "Pathotype" = metadata$Pathotype,
      "Resp" = metadata$DAS28.CRP.EULARresp.bin.V7,
      
      col = list(
        "CD3" = colorRamp2(c(0, 4), c("white", "slateblue")),
        "CD20"  = colorRamp2(c(0, 4), c("white", "gold2")),
        "CD68L"  = colorRamp2(c(0, 4), c("white", "darkcyan")),# hcl_palette = "Oranges", reverse = T),
        "CD68SL"  = colorRamp2(c(0, 4), c("white", "darkgreen")),
        "CD138"  = colorRamp2(c(0, 4), c("white", "darkorange")),
        "Pathotype" = c("Lymphoid" = 'blue',
                        "Myeloid" = 'red',
                        "Fibroid" = 'green', 
                        "Ungraded" = 'gray'),
        "Cell type" = c("Brich" = "darkblue", "Bpoor" = "lightblue"),
                   "Resp" = c("Good.or.mod" = "dodgerblue", 
                              "Non.responder" = "firebrick")
        ),
      annotation_name_side = "left",
      annotation_legend_param = list(direction = "horizontal", nrow = 1)
)

kmeans.rows <- kmeans(heatmap.df, centers = 3)
cluster_factor <- factor(kmeans.rows$cluster) 

hm <- Heatmap(heatmap.df, name = "z score", 
        top_annotation = ha,
        show_column_names = FALSE, show_row_names = FALSE,
        show_row_dend = FALSE,
        row_split = cluster_factor,
        cluster_row_slices = FALSE,
        column_split = 3,
        row_title = c("1", "2", "3"),  # to make things clearer
        clustering_distance_columns = "pearson",  #pearson
        clustering_method_columns = "complete", # ward.D2
        # column_title = "STRAP",
        heatmap_legend_param = list(direction = "vertical"),
        raster_by_magick = TRUE,
        raster_quality = 8,
        # use_raster = FALSE  # it'll make it slower, but necessary for very high quality 
)

# NOTE: if you run into troubles in trying to visualize this heatmap it might be  
# that you're having Cairo related issues. On a MAC, you might need to install
# XQuartz (https://www.xquartz.org/).

# if (plot.on.file) (svg("Plots/STRAP.heatmap.clusters.svg", height = 7, width = 5))
if (plot.on.file) (pdf("Plots/STRAP.heatmap.clusters_raster8.pdf", height = 7, width = 5))
# To make sure the order won't change:
hm1 <- draw(hm,
            show_heatmap_legend = FALSE,
            show_annotation_legend = FALSE) 
if (plot.on.file) dev.off()
# PS: the legend is plotted separately and is in common with the R4RA heatmap,
# see R4RA.clusters.R to plot the legend

# extract clusters of patients and 
# add annotation this annotation as new column of the metadata file
col_clusters <- column_order(hm1)
metadata[colnames(vst)[col_clusters[[1]]], "cluster"] <- "cluster 1"
metadata[colnames(vst)[col_clusters[[2]]], "cluster"] <- "cluster 2"
metadata[colnames(vst)[col_clusters[[3]]], "cluster"] <- "cluster 3"
saveRDS(metadata, "../Data/STRAP/STRAP.metadata_ML_clusters.RDS")

# extract and save gene clusters 
# (they will be used for pathway analysis and Venn diagrams)
row_clusters <- row_order(hm1)
names(row_clusters) <- c("1", "2", "3") # match with titles in plot
for (cluster in 1:3) {
  cluster.genes <- rownames(vst[row_clusters[[cluster]], ])
  write.csv(cluster.genes, paste0("STRAP.genes.cluster", cluster, ".csv"),
            row.names = F, quote = F)
}