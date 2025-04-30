rm(list = ls())
library(dplyr)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(circlize) 
library(amap)
library(Cairo)
plot.on.file <- F
ht_opt$message = FALSE # turn off message about use_raster

# load & filter data ---------------------------------------------------------------
metadata <- readRDS("../Data/R4RA/metadata.R4RA.RDS") %>% #baseline only
  filter(Visit == 3)
rownames(metadata) <- metadata$Seq_ID.V2
counts <- readRDS("../Data/R4RA/txi.R4RA.gene.symbols.RDS")$counts

coding.genes <- read.csv("../Data/gene_type.txt", stringsAsFactors = FALSE)[, c(1,2)] %>%
  filter(Gene.type == "protein_coding") %>%
  pull(Gene.name)
coding.genes <- coding.genes[coding.genes %in% rownames(counts)]

counts <- counts[which(rownames(counts) %in% coding.genes), rownames(metadata)]
keep <- filterByExpr(counts, min.count = 20, min.total.count = 300000)

table(keep)

counts <- counts[keep, ] # 2072 genes
nrow(counts)
mode(counts) <- "integer"

vst <- vst(counts, blind=TRUE) 

# Heatmap -----------------------------------------------------------------
# Rename levels for shorter legends 
levels(metadata$Cell.Type.V2) <- c("Bpoor", "Brich", "GC", "NA")
levels(metadata$DAS28.CRP.EULARresp.bin.V7) <- c("Good/Moderate", "No Resp")

heatmap.df <- t(scale(t(vst)))
heatmap.df <- heatmap.df[which(complete.cases(heatmap.df)),]

set.seed(10)

ha = HeatmapAnnotation(
                       "CD3" = metadata$CD3.V2,  
                       "CD20"  = metadata$CD20.V2,
                       "CD68L" = metadata$CD68L.V2,
                       "CD68SL" = metadata$CD68SL.V2,
                       "CD138"  = metadata$CD138.V2,
                       "Cell type" = metadata$Cell.Type.V2,
                       "Pathotype" = metadata$Pathotype.V2,
                       "Resp" = metadata$DAS28.CRP.EULARresp.bin.V7,
                       col = list(
                         "CD3" = colorRamp2(c(0, 4), c("white", "slateblue")),
                         "CD20"  = colorRamp2(c(0, 4), c("white", "gold2")),
                         "CD68L"  = colorRamp2(c(0, 4), c("white", "darkcyan")),
                         "CD68SL"  = colorRamp2(c(0, 4), c("white", "darkgreen")),
                         "CD138"  = colorRamp2(c(0, 4), c("white", "darkorange")),
                         "Cell type" = c("Brich" = "darkblue", "Bpoor" = "lightblue",
                                         "GC" = "darkviolet", "NA" = "gray"),
                         "Pathotype" = c("Lymphoid"='blue',
                                         "Myeloid"='red',
                                         "Fibroid"='green', 
                                         "Ungraded" = 'gray'),
                         "Resp" = c("Good/Moderate" = "dodgerblue",
                                             "No Resp" = "firebrick")
                         ),
                       annotation_name_side = "left",
                       annotation_legend_param = list(direction = "vertical",
                                                      nrow = 1)
)

col_dend = as.dendrogram(hclust(Dist(t(heatmap.df), method = "pearson"),
                                method = "complete"))

# generate and reorder column dendogram to have the fibroid cluster in the middle
# (for coherence with STRAP)
tmp <- col_dend[[2]][[1]]
col_dend[[2]][[1]] <- col_dend[[2]][[2]]
col_dend[[2]][[2]] <- tmp

# generate and reorder row clusters 
# (for coherence with STRAP)
kmeans.rows <- kmeans(heatmap.df, centers = 3)
cluster_factor <- factor(kmeans.rows$cluster, levels = c(2, 3, 1))

hm <-
  Heatmap(heatmap.df, name = "z score", 
          top_annotation = ha,
          show_column_names = FALSE, show_row_names = FALSE,
          show_row_dend = FALSE,
          row_split = cluster_factor,
          cluster_columns = col_dend,
          column_split = 3,
          row_title = c("1", "2", "3"),  # to make things clearer
          heatmap_legend_param = list(direction = "horizontal"),
          raster_by_magick = TRUE,
          raster_quality = 2
          # use_raster = FALSE, # for maximum quality of the image
)

# if (plot.on.file) (svg("Plots/R4RA.heatmap.clusters.svg", height = 7, width = 5))
if (plot.on.file) (pdf("Plots/R4RA.heatmap.clusters_raster2.pdf", height = 7, width = 5))
# To make sure the order won't change 
hm1 <- draw(hm, annotation_legend_side = "bottom",
            show_heatmap_legend = FALSE,
            show_annotation_legend = FALSE
            ) 
if (plot.on.file) dev.off()

# plot legend as separate figure (in common with STRAP)
hm2 <- draw(hm, annotation_legend_side = "right",
            show_heatmap_legend = TRUE,
            show_annotation_legend = TRUE)

# if (plot.on.file) (svg("Plots/R4RA.anno.heatmap.legend.svg", height = 8, width = 9.5))
if (plot.on.file) (pdf("Plots/R4RA.anno.heatmap.legend.pdf", height = 8, width = 9.5))
grid.newpage()
draw_annotation_legend(hm2)
draw_heatmap_legend(hm2)
if (plot.on.file) dev.off()


# extract clusters of patients and 
# add annotation this annotation as new column of the metadata file
col_clusters <- column_order(hm1)
col_clusters <- list(col_clusters[[1]], col_clusters[[3]], col_clusters[[2]])
col_cluster_factor <- factor(kmeans.rows$cluster, levels = c(2, 3, 1))

metadata[colnames(vst)[col_clusters[[1]]], "cluster"] <- "cluster 1"
metadata[colnames(vst)[col_clusters[[2]]], "cluster"] <- "cluster 2"
metadata[colnames(vst)[col_clusters[[3]]], "cluster"] <- "cluster 3"
saveRDS(metadata, "metadata.R4RA_clusters.RDS")

# extract and save gene clusters 
# (they will be used for the Venn diagrams)
row_clusters <- row_order(hm1)
names(row_clusters) <- c("1", "2", "3") # match with titles in plot
for (cluster in 1:3) {
  cluster.genes <- rownames(vst[row_clusters[[cluster]], ])
  write.csv(cluster.genes, paste0("R4RA.genes.cluster", cluster, ".csv"),
            row.names = F, quote = F)
}
