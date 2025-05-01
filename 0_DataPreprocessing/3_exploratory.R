rm(list=ls())

library(patchwork)
library(ggplot2)

#### QC checks for outliers

legendtitle <- ""
load("STRAP_db_v202103.RData")
metadata <- STRAP_db_v202103$RNAseq_samples
vst <- STRAP_db_v202103$RNAseq_data$vst

scale <- T
norm.vals <- vst
dim(norm.vals)
norm.vals <- norm.vals[which(!apply(norm.vals,1, function(x) all(x<=2))),]
dim(norm.vals)
norm.vals <- norm.vals[which(!apply(norm.vals,1, function(x) all(x==x[1]))),]
dim(norm.vals)
pc <- prcomp(t(norm.vals), scale. = scale)
varExp <- round(pc$sdev^2/sum(pc$sdev^2),digits = 2) * 100
summary(pc)


pc_df <- as.data.frame(pc$x)
pc_df$group <- as.numeric(metadata$MappingRatesSalmon)
pc_df$Outlier <- metadata$Outliers
p<-ggplot(pc_df,aes(x=PC1,y=PC2))
p <- p + geom_point(aes(shape=Outlier, color=group), size=4, alpha = 0.7)
p <- p + theme_classic2() + scale_color_gradient(low = "#f8e0d9", high = "#c34924", na.value = "grey") 
p <- p + xlab(paste0("PC1 variance %", varExp[1]))
p <- p + ylab(paste0("PC2 variance %", varExp[2]))
p <- p + labs(col=legendtitle) + theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(colour="Mapping rates", title = ifelse(scale, "All samples\n", "All samples\n")) +
    theme(legend.position="bottom", legend.box = "horizontal",  legend.direction = "vertical", legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "gray96", color = NA), legend.justification = c(0, 1)) +
  geom_vline(aes(xintercept=-155), color="grey50", linetype="dashed", size=1) + 
    theme(axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"))

p
