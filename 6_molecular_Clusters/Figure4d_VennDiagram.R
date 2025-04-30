rm(list = ls())
library(grid)
library(VennDiagram)
library(polyclip)

plot.on.file <- F
cluster <- 1  # 1, 2, or 3

STRAP.genes <- read.csv(paste0("STRAP.genes.cluster", cluster, ".csv"))[,1]
R4RA.genes  <- read.csv(paste0("R4RA.genes.cluster", cluster, ".csv"))[,1]

grid.newpage()
venn <- venn.diagram(x = list(R4RA.genes, STRAP.genes), alpha = 0.8,
                     category.names = c("STRAP", "R4RA"),
                     main = paste0("Cluster ", cluster),
                     cat.fontfamily = "sans", cat.cex = 3, cat.dist = c(0.17,0.17),
                     fontfamily = "sans", cex = 3,
                     filename = NULL,
                     print.mode = "raw",  
                     fill = c("gold3", "blue"),
                     lwd = 2,
                     disable.logging = T,
                     scaled = FALSE
)

# To customise the colors of the two circle and their intersection I need some 
# to plot each element separately. 

# Extract the coordinates for the R4RA and STRAP sets from the venn object
R4RA  <- list(list(x = as.vector(venn[[3]][["x"]]), y = as.vector(venn[[3]][["y"]])))
STRAP <- list(list(x = as.vector(venn[[4]][["x"]]), y = as.vector(venn[[4]][["y"]])))

# Calculate the intersection between R4RA and STRAP sets using polyclip
intersection <- polyclip(R4RA, STRAP)

# Extract the text labels from the venn object
# First identify which elements in venn contain 'text' in their names
ix <- sapply(venn, function(x) grepl("text", x$name, fixed = TRUE))
# Then extract the coordinates and labels for these text elements
labs <- do.call(rbind.data.frame, lapply(venn[ix], `[`, c("x", "y", "label")))

# if(plot.on.file)(svg(filename = paste0("Plots/Venn.cluster", cluster, ".svg"),
#                                       width = 4, height = 4))
if(plot.on.file)(pdf(paste0("Plots/Venn.cluster", cluster, ".pdf"),
                     width = 4, height = 4))

# Create an empty plot with no axes
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
# Draw the R4RA set as a blue polygon
polygon(R4RA[[1]], border = NA, col = adjustcolor("blue", alpha.f=0.8))
# Draw the STRAP set as a gold polygon
polygon(STRAP[[1]], border = NA, col = adjustcolor("gold3", alpha.f=0.8))
# Draw the intersection as a green polygon
polygon(intersection[[1]], border = NA, col = adjustcolor("green3", alpha.f=0.5))
# Add text labels to the diagram
text(x = labs$x[1:3], y = labs$y[1:3], labels = labs$label[1:3], cex = 1.8, col = "white")
text(x = labs$x[4:5], y = labs$y[4:5], labels = labs$label[4:5], cex = 1.5)
text(x = labs$x[6], y = labs$y[6], labels = labs$label[6], cex = 1.5)

if(plot.on.file)(dev.off())
