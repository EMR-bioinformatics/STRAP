# NOTE: all DESeq2 DEGs files come from figure 1

rm(list = ls())
library(volcano3D)
library(DESeq2)
library(plotly)
library(devtools)
library(dplyr)
library(Cairo)

############### generate polar objects ##################
res_eta <- readRDS("deseq2.DEGs_Etanercept.rds")
res_rtx <- readRDS("deseq2.DEGs_Rituximab.rds")
res_toc <- readRDS("deseq2.DEGs_Tocilizumab.rds")

rn <- unique(c(rownames(res_eta), rownames(res_rtx), rownames(res_toc)))

df1 <- data.frame(Eta = res_eta[rn, 'log2FoldChange'], # invert sign here to switch between resp/non resp. For non-responders: -res_eta[rn, 'log2FoldChange'],
                  Rtx = res_rtx[rn, 'log2FoldChange'], # invert sign here to switch between resp/non resp
                  Toc = res_toc[rn, 'log2FoldChange'], # invert sign here to switch between resp/non resp
                  SE_Eta = res_eta[rn, 'lfcSE'],
                  SE_Rtx = res_rtx[rn, 'lfcSE'],
                  SE_Toc = res_toc[rn, 'lfcSE'], row.names = rn)

df1 <- volcano3D:::polar_xy(df1)

pvals <- data.frame(Eta = res_eta[rn, 'pvalue'],
                    Rtx = res_rtx[rn, 'pvalue'],
                    Toc = res_toc[rn, 'pvalue'], row.names = rn)
pvals <- as.matrix(pvals)

padj <- data.frame(Eta = res_eta[rn, 'qvalue'],
                   Rtx = res_rtx[rn, 'qvalue'],
                   Toc = res_toc[rn, 'qvalue'], row.names = rn)
padj <- as.matrix(padj)

scheme <- c('grey60', 'red', 'gold2', 'green3', 
            'cyan', 'blue', 'purple', 'black')
outcome <- factor(levels=(c("Eta", "Rtx", "Toc")))


polar_p2 <- function(df, pvals, padj = pvals, pcutoff = 0.05,
                     scheme = c('grey60', 'red', 'gold2', 'green3', 
                                'cyan', 'blue', 'purple', 'black'),
                     labs = NULL) {
  outcome_levels <- colnames(df)[1:3]
  pcheck <- sapply(1:3, function(i) padj[,i] < pcutoff & df[,i] > 0)
  colnames(pcheck) <- outcome_levels
  pcheck <- pcheck *1  # convert to numeric
  pcheck[is.na(pcheck)] <- 0
  # negatives
  # pneg <- sapply(1:3, function(i) padj[,i] < pcutoff & df[,i] < 0)
  # colnames(pneg) <- outcome_levels
  # negrs <- rowSums(pneg, na.rm=T)
  # pneg[negrs>0] <- 1 - pneg[negrs>0]  # invert
  # pneg[is.na(pneg)] <- 0
  # posrs <- rowSums(pcheck, na.rm=T)
  # pcheck[posrs == 0, ] <- pneg[posrs == 0, ]  # replace only grey points with pneg
  pset <- apply(pcheck, 1, paste, collapse = "")  # slow
  
  if (is.null(labs) | length(labs) == 3) {
    abbrev <- if (length(labs) == 3) labs else abbreviate(outcome_levels, 1)
    labs <- c("ns",
              paste0(abbrev[1], "+"),
              paste0(abbrev[1], "+", abbrev[2], "+"),
              paste0(abbrev[2], "+"),
              paste0(abbrev[2], "+", abbrev[3], "+"),
              paste0(abbrev[3], "+"),
              paste0(abbrev[1], "+", abbrev[3], "+"),
              paste0(abbrev[1], "+", abbrev[2], "+", abbrev[3], "+"))
  }
  lab <- factor(pset, levels = c('000', '100', '110', '010',
                                 '011', '001', '101', '111'),
                labels = labs)
  col <- scheme[lab]
  data.frame(col, lab)
}

ptab <-polar_p2(df1, pvals, padj)

table(ptab$lab, useNA = "always")

df1$col <- df1$lab <- NULL
df1 <- cbind(df1, ptab)

custom_polar <- methods::new("volc3d",
                             df = list(scaled = df1, unscaled = NULL),
                             outcome = outcome,
                             data = data.frame(), pvals = pvals, padj = padj,
                             pcutoff = 0.05, scheme = scheme,
                             labs = levels(ptab$lab))

saveRDS(custom_polar, "polar.object.DEGs.from.Fig1_Resp.RDS")
# saveRDS(custom_polar, "polar.object.DEGs.from.Fig1_nonResp.RDS")

############# plot ################
# figure 2c and 2d
custom_polar <- readRDS("polar.object.DEGs.from.Fig1_Resp.RDS")

sig.genes <- rownames(custom_polar@df$scaled[which(custom_polar@df$scaled$lab != "ns"), ])

radial_plotly(custom_polar,
              arrow_length = 4,
              axis_angle = 1/6) %>%
  config(edits = list(annotationTail = TRUE),
         toImageButtonOptions = list(format = "svg")) %>%
  # toWebGL() %>%   # add webGL for faster, but lower quality figures 
  plotly::layout(autosize = F, width = 500, height = 500)
