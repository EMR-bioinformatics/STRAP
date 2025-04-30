library(caret)
library(glmnet)
library(gbm)
library(xgboost)
library(ggplot2)
library(dplyr)
library(randomForest)
library(mda)
library(pls)
library(MLmetrics)
library(RhpcBLASctl)
library(devtools)
library(CORElearn)
library(nestedcv) #package 0.7.9
library(patchwork)
library(ggpubr)

library(parallel)

meta_BL <- readRDS("R4RA_STRAP_meta.rds")

RhpcBLASctl::omp_set_num_threads(1L) #switch off multithreading for xgboost to work 
source("compare_mod.R")

#curation of vst on R4RA + STRAP rawcounts

load("merge_R4RA_STRAP.RData")

all(colnames(rawcounts_all) == meta_RNA_samples$SeqID)
#TRUE

rownames(meta_RNA_samples) <- meta_RNA_samples$SeqID

library(DESeq2)

dds <- DESeqDataSetFromMatrix(round(rawcounts_all),
                              meta_RNA_samples,
                              design= formula(~Cohort))

vsd <- vst(dds, blind = FALSE)
vst_all <- assay(vsd)

vst_BL <- vst_all[ , colnames(vst_all) %in% meta_BL$SeqID]

#all(colnames(vst_BL) == meta_BL$SeqID)
#TRUE

#subset to nanostring 
nanostring <- readRDS("geneAnnot_nCounter_custompanel2.rds")

nanostring <- nanostring[nanostring$V1 == "Endogenous", ]
n_vst <- intersect(nanostring$V2, rownames(vst_BL))

data <- t(vst_BL[n_vst, ]) 

#colnames(data)[nearZeroVar(data)]
#"IFNA8" "IL9"   "IFNA2"

data2 <- data[, -nearZeroVar(data)]

rtx <- meta_BL$Treatment == "Rituximab"
toc <- meta_BL$Treatment == "Tocilizumab"

# original
oy <- as.vector(meta_BL$TargetDAS28.ESR.V7)
oy <- factor(make.names(oy))
names(oy) <- meta_BL$SeqID

clinset2 <- c("sqTJ", "sqSJ", "sqESR", "lnCRP", "Arthritis.Activity", "CCP", "RF")
#rename for ease/later plotting
colnames(meta_BL)[colnames(meta_BL) == "logCRP"] <- "lnCRP"
rownames(meta_BL) <- meta_BL$SeqID
clindata2 <- meta_BL[ , clinset2]

#restrict nanostring data
nanostring_mean <- colMeans(data2)
lowgene <- names(nanostring_mean)[nanostring_mean < 6]
data3 <- data2[ , !colnames(data2) %in% lowgene]

data3clin <- cbind(data3, clindata2)

####nestedcv####

#rituximab

grd <- list(
  filterFUN = "ttest_filter",
  filterFUN = "glmnet_filter",
  filterFUN = "pls_filter")

trc <- trainControl(method = "cv",
                    number = 10,
                    classProbs = TRUE,
                    savePredictions = "final",
                    summaryFunction = defaultSummary)

mods <- list(
  glmnet = list(model="glmnet", args=list(family="binomial",
                                          alphaSet = seq(0.7, 1, 0.1),
                                          filter_options = list(nfilter = 40),
                                          min_1se = 0)),
  rf = list(model="rf", args=list(filter_options = list(nfilter = 40),
                                  balance="smote")),
  gbm = list(model="gbm", args=list(metric="Kappa",
                                    trControl = trc,
                                    filter_options = list(nfilter = 40))),
  mda = list(model="mda", args=list(metric="Kappa",
                                    trControl = trc,
                                    filter_options = list(nfilter = 40))),
  pls = list(model = "pls", args=list(metric="Kappa",
                                      trControl = trc,
                                      filter_options = list(nfilter = 40))),
  svmRadial = list(model="svmRadial", args=list(metric="Kappa",
                                                trControl = trc,
                                                filter_options = list(nfilter = 40),
                                                balance="smote")),
  svmPoly = list(model="svmPoly", args=list(metric="Kappa",
                                            trControl = trc,
                                            filter_options = list(nfilter = 40),
                                            balance="smote")),
  xgbTree = list(model="xgbTree", args=list(metric="Kappa",
                                            trControl = trc,
                                            filter_options = list(nfilter = 40))),
  xgbLinear = list(model="xgbLinear", args=list(metric="Kappa",
                                                trControl = trc,
                                                filter_options = list(nfilter = 40)))
)

resRtx <- compare_mod(oy[rtx], data3clin[rtx, ], mods, grid = grd,
                      repeats = 25, n_outer_folds = 10, rep.cores = 15,
                      return_fits = T)

# saveRDS(resRtx, "Final/resRtx.rds")


resToc <- compare_mod(oy[toc], data3clin[toc, ], mods, grid = grd,
                      repeats = 25, n_outer_folds = 10, rep.cores = 13,
                      return_fits = T)

#saveRDS(resToc, "Final/resToc.rds")

pp <- wrap_plots(
  plot_modcompare(resToc, xvar2 = "filter") + ggtitle("Tocilizumab"),
  plot_modcompare(resRtx, xvar2 = "filter") + ggtitle("Rituximab"),
  guides = "collect") &
  labs(x = "") 

ggsave("compare_mods_AUC.pdf", pp, width = 9, height = 3.5)

####nfilter####

nfilter_grd <- list(filter_options = list(nfilter = 15),
                    filter_options = list(nfilter = 20),
                    filter_options = list(nfilter = 30),
                    filter_options = list(nfilter = 40),
                    filter_options = list(nfilter = 50))

pls_mod <- list(pls = list(model = "pls", args=list(filterFUN = "ttest_filter",
                                                    metric="Kappa",
                                                    trControl = trc)))

resRtx_tune <- compare_mod(oy[rtx], data3clin[rtx, ], pls_mod, grid = nfilter_grd,
                      repeats = 25, n_outer_folds = 10, rep.cores = 15,
                      return_fits = T)

# saveRDS(resRtx_tune, "resRtx_tune.rds")

glmnet_mod <- list(glmnet = list(model="glmnet", args=list(filterFUN = "ttest_filter",
                                                           family="binomial",
                                                           alphaSet = seq(0.7, 1, 0.1),
                                                           min_1se = 0)))

resToc_tune <- compare_mod(oy[toc], data3clin[toc, ], glmnet_mod, grid = nfilter_grd,
                           repeats = 25, n_outer_folds = 10, rep.cores = 15,
                           return_fits = T)

# saveRDS(resToc_tune, "resToc_tune.rds")

pp3 <- wrap_plots(
  plot_modcompare(resToc_tune, xvar2 = "nfilter") + ggtitle("Tocilizumab"),
  plot_modcompare(resRtx_tune, xvar2 = "nfilter") + ggtitle("Rituximab"),
  guides = "collect") &
  labs(x = "") 

ggsave("nfilter_AUC.pdf", pp3, width = 9, height = 3.5)

####set seed model####

set.seed(11, "L'Ecuyer-CMRG")

rtxMod <- nestcv.train(oy[rtx], data3clin[rtx, ],
                       filterFUN = "ttest_filter",
                       method = "pls", 
                       metric="Kappa",
                       trControl = trc,
                       filter_options = list(nfilter = 20),
                       cv.cores = 8)

summary(rtxMod)
#AUC 0.7499 BA 0.6919 

#saveRDS(rtxMod, "rtxMod.rds")

set.seed(23, "L'Ecuyer-CMRG")

tocMod <- nestcv.glmnet(oy[toc], data3clin[toc, ],
                        filterFUN = "ttest_filter",
                        family="binomial",
                        alphaSet = seq(0.7, 1, 0.1),
                        min_1se = 0,
                        filter_options = list(nfilter = 30),
                        cv.cores = 8)

summary(tocMod)
#AUC  0.7850 BA 0.7125 

#saveRDS(tocMod, "tocMod.rds")

#var stability/importance
v1 <- plot_var_stability(tocMod, direction = 1, scheme = c("red", "royalblue")) +
  ggtitle("Tocilizumab")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

v2 <- plot_var_stability(rtxMod, direction = 1, scheme = c("red", "royalblue"),
                         dir_labels = c("Up in Non-responder", "Up in Responder")) +
  ggtitle("Rituximab")+
  theme(legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill="white",
                                         colour="white"),
        plot.title = element_text(hjust = 0.5))

pp4 <- wrap_plots(v1, v2)

ggsave("varImp.pdf", pp4, width = 8, height = 7)

#roc curves

r1 <- innercv_roc(tocMod)
r2 <- innercv_roc(rtxMod)

pdf("final_roc.pdf",
    width=9.5, height=3.2)
op <- par(pty="s", mfrow=c(1, 3), tcl=-0.4)
plot(r1, col="gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Tocilizumab model (glmnet)", font.main = 1)
lines(tocMod$roc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(tocMod$roc$auc, r1$auc), digits=3)),
       bty="n")
plot(r2, col="gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Rituximab (pls)", font.main = 1)
lines(rtxMod$roc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(rtxMod$roc$auc, r2$auc), digits=3)),
       bty="n")
par(op)
dev.off()

####plot variable ranks plot####

source("var_ranks.R")

pp4 <- wrap_plots(plot_var_ranks(tocMod) + ggtitle("Tocilizumab"), 
                  plot_var_ranks(rtxMod) + ggtitle("Rituximab")) &
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Final/var_ranks.pdf", pp4, width = 8, height = 7)
