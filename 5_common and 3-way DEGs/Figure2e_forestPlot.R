rm(list = ls())
library(dplyr)
library(volcano3D)
library(gtools)
library(ggplot2)
library(Cairo)

resp.polar <- readRDS("polar.object.DEGs.from.Fig1_Resp.RDS")
# nonResp.polar <- readRDS("polar.object.DEGs.from.Fig1_nonResp.RDS")

# let's pick some example:
genes <- c("KRT10", "MMP9", "MS4A1", "PAX5", "CR2")

# Figure 2d (forest plot) -------------------------------------------------
forest.df <- resp.polar@df[[1]][genes,1:3] # using nonResp.polar would only change the direction
forest.df$id <- rownames(forest.df)
CI.df <- resp.polar@df[[1]][genes,4:6] * 1.96 # I want to show 95% confidence intervals 
padj.df  <- as.data.frame(resp.polar@padj)[genes, ]

forest.df <- tidyr::gather(forest.df, key = "drug", value = "x", 1:3)
CI.df     <- tidyr::gather(CI.df,     key = "drug", value = "CI", 1:3)
padj.df   <- tidyr::gather(padj.df,   key = "drug", value = "padj", 1:3)

forest.df$CI <- CI.df$CI
forest.df$LCI <- forest.df$x - forest.df$CI
forest.df$UCI <- forest.df$x + forest.df$CI
forest.df$padj <- padj.df$padj
forest.df$padj.star <- stars.pval(forest.df$padj)
forest.df$drug <- as.factor(forest.df$drug)
forest.df$gene <- factor(forest.df$id, levels = genes)
forest.df$id <- forest.df$drug
forest.df$drug <- NULL

vdiff <- diff(range(c(forest.df$LCI, forest.df$UCI, 0), na.rm=TRUE))

data <- forest.df
xmax <- max(data$x + data$CI, na.rm = TRUE)
write.csv(data, "forest.table.csv")

ggplot(data, aes(y = .data$id, x = .data$x, color = .data$id)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = .data$x - .data$CI, xmax = .data$x + .data$CI),
                width = 0.3) + 
  geom_point() +
  scale_color_manual(values = c("red", "green3", "blue")) + 
  facet_grid(gene ~ .) +
  xlab(bquote("log"[2] ~ "FC")) + ylab("") + labs(color = "") +
  geom_text(y = data$id, x = xmax + vdiff * 0.08, colour = "black", 
            label = data$padj.star, show.legend = FALSE) +
  expand_limits(x = xmax + vdiff * 0.1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "cm"), 
        strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0)) + 
  coord_cartesian(clip = 'off') +
  annotation_custom(grid::textGrob(x = unit(0.8, 'npc'),
                                   y = unit(-6, 'mm'),
                                   label = "Responders",
                                   gp = grid::gpar(col = 'black', cex = 0.6,
                                                   fontfamily = 'sans')) ) +
  annotation_custom(grid::linesGrob(x = unit(c(0.68, 0.94), 'npc'),
                                    y = unit(c(-8, -8), 'mm'),
                                    arrow = arrow(length = unit(1.5, 'mm')),
                                    gp = grid::gpar(col = 'black')) ) +
  annotation_custom(grid::textGrob(x = unit(0.2, 'npc'),
                                   y = unit(-6, 'mm'),
                                   label = "Non Responders",
                                   gp = grid::gpar(col = 'black', cex = 0.6,
                                                   fontfamily = 'sans')) ) +
  annotation_custom(grid::linesGrob(x = unit(c(0.35, 0.02), 'npc'),
                                    y = unit(c(-8, -8), 'mm'),
                                    arrow = arrow(length = unit(1.5, 'mm')),
                                    gp = grid::gpar(col = 'black')) )
ggsave("genes.forest.pdf", device = "pdf", width = 3.6, height = 3.3) 
