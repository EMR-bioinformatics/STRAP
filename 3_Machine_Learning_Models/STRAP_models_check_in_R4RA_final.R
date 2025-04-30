library(nestedcv)
library(caret)
library(gbm)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(SuperLearner)
library(pROC)
source("/users/myles/documents/github/ML_performance/compare_mod.R")

load("/Users/myles/R/STRAP/STRAP_more_models_290424.RData")

nanostring <- readRDS("/Users/myles/R/STRAP/geneAnnot_nCounter_custompanel2.rds")
nanostring <- nanostring[nanostring$V1=="Endogenous", ]
n_vst <- intersect(nanostring$V2, rownames(vstdata))

data <- t(vstdata[n_vst, ])

eta <- visit3 & metadata$Randomised.Medication == "Etanercept"
rtx <- visit3 & metadata$Randomised.Medication == "Rituximab"
toc <- visit3 & metadata$Randomised.Medication == "Tocilizumab"

metadata$sqTJ <- sqrt(metadata$Number.of.Tender.Joints)
metadata$sqSJ <- sqrt(metadata$Number.of.Swollen.Joints)
metadata$sqESR <- sqrt(metadata$ESR)
metadata$lnCRP <- log(metadata$CRP +1)

data2 <- data

# caret function for finding near zero vars
# colnames(data2)[nearZeroVar(data2)]
data2 <- data2[, -nearZeroVar(data2)]
dim(data2)

# original
oy <- metadata$TargetDAS28.ESR.V7
oy <- factor(make.names(oy))
names(oy) <- rownames(metadata)

# use 4 level DAS28.ESR.status
y <- metadata$DAS28.ESR.status.V7
y <- (as.numeric(y) - 1) / 3
# hist(y[rtx])

# use continuous DAS28.ESR.V7
out_v7 <- Outcome[Outcome$Visit == 7, ]
metadata$DAS28.ESR.V7 <- out_v7$DAS28.ESR[match(metadata$Patient.I.D., out_v7$Patient.I.D.)]
# plot(y = metadata$DAS28.ESR.V7, x = metadata$DAS28.ESR.status.V7)
# hist(-metadata$DAS28.ESR.V7)
# y <- metadata$DAS28.ESR.V7

# expanded geneset for RTX
clinset2 <- c("sqTJ", "sqSJ", "sqESR", "lnCRP", "CCP", "RF", "Arthritis.Activity")
clinset3 <- c("sqTJ", "sqSJ", "sqESR", "lnCRP", "CCP", "RF")
clindata2 <- metadata[, clinset2]
clindata3 <- metadata[, clinset3]

# load final models
load("/Users/myles/Dropbox/R scripts/STRAP/STRAP models/STRAP_final_models.Rdata")

# load merged dataset
load("/Users/myles/R/merged/merge_R4RA_STRAP.RData")
dim(meta_RNA_samples)
dim(vst_all_comBat)

# check numbers
ind0 <- meta_RNA_samples$Cohort == "R4RA" & meta_RNA_samples$Visit == 3
table(meta_RNA_samples$Treatment[ind0])

# check tocMod (STRAP) in R4RA
all(tocMod$final_vars %in% rownames(vst_all_comBat))
tocMod$final_vars[!tocMod$final_vars %in% rownames(vst_all_comBat)]
ind <- meta_RNA_samples$Cohort == "R4RA" & meta_RNA_samples$Visit == 3 &
  meta_RNA_samples$Treatment == "Tocilizumab"
vst_all_comBat_clin <- t(vst_all_comBat)
clin <- data.frame(sqTJ = sqrt(meta_RNA_samples$Tender),
                   sqSJ = sqrt(meta_RNA_samples$Swollen),
                   Arthritis.Activity = meta_RNA_samples$Pain)
vst_all_comBat_clin <- cbind(vst_all_comBat_clin, clin)

r4ra_toc_dat <- vst_all_comBat_clin[ind, tocMod$final_vars]

r4ra_toc_preds <- predict(tocMod, r4ra_toc_dat)
true_y <- meta_RNA_samples$TargetDAS28.ESR.V7[ind]
true_y <- factor(true_y, labels = c("Non.responder", "Responder"))
table(true_y, r4ra_toc_preds)
confusionMatrix(table(true_y, r4ra_toc_preds), positive = "Responder")

r4ra_toc_preds <- predict(tocMod, r4ra_toc_dat, type = "prob")
r4ra_toc_roc <- pROC::roc(true_y, r4ra_toc_preds[,2], direction = "<")
pROC::auc(r4ra_toc_roc)

r4ra_toc_youden <- pROC::coords(r4ra_toc_roc, x = "best")$threshold
pred_y <- r4ra_toc_preds[,2] > r4ra_toc_youden
pred_y <- factor(pred_y, labels = c("Non.responder", "Responder"))

table(true_y, pred_y)
cm <- confusionMatrix(table(true_y, pred_y), positive = "Responder")

# check rtxMod (STRAP) in R4RA
ind <- meta_RNA_samples$Cohort == "R4RA" & meta_RNA_samples$Visit == 3 &
  meta_RNA_samples$Treatment == "Rituximab"
vst_all_comBat_clin <- t(vst_all_comBat)
clin <- data.frame(sqTJ = sqrt(meta_RNA_samples$Tender),
                   sqESR = sqrt(meta_RNA_samples$ESR),
                   Arthritis.Activity = meta_RNA_samples$Pain)
vst_all_comBat_clin <- cbind(vst_all_comBat_clin, clin)

r4ra_rtx_dat <- vst_all_comBat_clin[ind, rtxMod$final_vars]
r4ra_rtx_preds <- predict(rtxMod, r4ra_rtx_dat)
true_y <- meta_RNA_samples$TargetDAS28.ESR.V7[ind]

r4ra_rtx_roc <- pROC::roc(true_y, r4ra_rtx_preds, direction = "<")
pROC::auc(r4ra_rtx_roc)
r4ra_rtx_youden <- pROC::coords(r4ra_rtx_roc, x = "best")$threshold

pred_y <- r4ra_rtx_preds > 0.5
pred_y <- factor(pred_y, labels = c("Non.responder", "Responder"))
true_y <- factor(true_y, labels = c("Non.responder", "Responder"))

table(true_y, pred_y)
confusionMatrix(table(true_y, pred_y), positive = "Responder")

# plot ROC curves (STRAP models tested in R4RA)
pdf("check_mods_R4RA_v2b.pdf",
    width = 8, height = 4)
op <- par(mfrow = c(1, 2), tcl = -0.3)
plot(r4ra_toc_roc, las = 1, mgp = c(2, 0.5, 0), col = "red",
     main = "Tocilizumab model tested in R4RA", font.main = 1)
legend("bottomright", legend = paste("AUC:", signif(pROC::auc(r4ra_toc_roc), 3)),
       bty = "n")
plot(r4ra_rtx_roc, las = 1, mgp = c(2, 0.5, 0), col = "#00abab",
     main = "Rituximab model tested in R4RA", font.main = 1)
legend("bottomright", legend = paste("AUC:", signif(pROC::auc(r4ra_rtx_roc), 3)),
       bty = "n")
par(op)
dev.off()

# replot STRAP ROC curves
plot(innercv_roc(tocMod), col = "gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Tocilizumab model (GBM)", font.main = 1)
lines(tocMod$roc, col = "blue")

##########
# original model building code
correl_filter(y[rtx], clindata2[rtx,], type="full", method = "spearman")
ttest_filter(y[rtx], clindata2[rtx,], type="full")

nz <- nearZeroVar(t(vstdata))
data3 <- t(vstdata[-nz, ])
cm <- matrixStats::colMedians(data3)
# cm <- matrixStats::colQuantiles(data3[rtx, ], probs = 0.75)
# hist(cm)
table(cm > 8)
data3 <- data3[, cm > 8]  # 3 gives good result but lots of weird genes
data3 <- cbind(data3, clindata2)
data3 <- as.matrix(data3)

ind <- removeAC(data3)
sum(ind)
# colnames(data3)[ind]
data3 <- data3[, !ind]
dim(data3)

# clinical only model
rtxClin <- nestcv.train(y[rtx], clindata2[rtx, ],
                       method = "rf",
                       cv.cores = 10)
summary(rtxClin)
plot_var_stability(rtxClin, direction = 2)

# RTX model based on regression
data2clin <- cbind(data2, clindata2)
dim(data2clin) # 273 503

set.seed(8, "L'Ecuyer-CMRG")
rtxMod <- nestcv.train(y[rtx], data2clin[rtx, ],
                       method = "xgbLinear",
                       filterFUN = glmnet_filter,
                       filter_options = list(nfilter = 28,
                                             force_vars = c("sqTJ", "sqESR", "Arthritis.Activity")),
                       cv.cores = 10)
summary(rtxMod)  # RMSE 0.3054  Rsq 0.2409

rtxMod$output$oy <- oy[match(rownames(rtxMod$output), names(oy))]
rtxroc <- roc(rtxMod$output$oy, rtxMod$output$predy, direction="<")
incv <- innercv_preds(rtxMod)
incv$oy <- factor(incv$testy > 0.5, labels = c("Non.responder", "Responder"))
rtx.inroc <- roc(incv$oy, incv$predy, direction = "<")
plot(rtxroc, col = "royalblue", las = 1)
lines(rtx.inroc, col = "gold")
rtxroc$auc  # 0.7535
rtx.inroc$auc  # 0.7335

plot_var_stability(rtxMod, direction = 2)
plot_var_stability(rtxMod, direction = 1, dir_labels = c("Up in Non.responder", "Up in Responder"))
plot(rtxMod$output$testy, rtxMod$output$predy)

linmods <- list(
  glmnet = list(model="glmnet", args=list(family="gaussian",
                                          filterFUN = "glmnet_filter",
                                          filter_options = list(nfilter = 30),
                                          alphaSet = seq(0.8, 1, 0.1),
                                          min_1se = 0)),
  rf = list(model="rf", args=list(filterFUN = "glmnet_filter",
                                  filter_options = list(nfilter = 30))),
  gbm = list(model="gbm", args=list(filterFUN = "glmnet_filter",
                                    filter_options = list(nfilter = 30))),
  # mda = list(model="mda", args=list(filterFUN = "lm_filter",
  #                                   filter_options = list(nfilter = 40))),
  svmRadial = list(model="svmRadial", args=list(filterFUN = "glmnet_filter",
                                                filter_options = list(nfilter = 30))),
  svmPoly = list(model="svmPoly", args=list(filterFUN = "glmnet_filter",
                                            filter_options = list(nfilter = 30))),
  xgbTree = list(model="xgbTree", args=list(filterFUN = "glmnet_filter",
                                            filter_options = list(nfilter = 30))),
  xgbLinear = list(model="xgbLinear", args=list(filterFUN = "glmnet_filter",
                                                filter_options = list(nfilter = 30))),
  SuperLearner = list(model="SL",
                      args=list(SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
                                filterFUN = "glmnet_filter",
                                filter_options = list(nfilter = 30)))
)
# mda sometimes fails to fit (errors)
# also try SuperLearner

resRtx <- compare_mod(y[rtx], data2clin[rtx, ], linmods,
                      repeats = 10, n_outer_folds = 5, rep.cores = 10)
resRtx$output
summary(resRtx)
plot_modcompare(resRtx, "RMSE")
plot_modcompare(resRtx, "Rsquared")

saveRDS(resRtx, "/../STRAP/resRtx.rds")

plot_modcompare(resRtx, "AUC")


##########################
# rebuild Eta, Toc models
data2 <- data2[, colnames(data2) != "ERCC3"]
y <- metadata$TargetDAS28.ESR.V7
y <- factor(make.names(y))

# ETA
set.seed(12, "L'Ecuyer-CMRG")  # AUC 0.7647, BA 0.731
etaMod <- nestcv.glmnet(y[eta], data2[eta, ],
                        family="binomial",
                        filterFUN = ttest_filter,
                        filter_options = list(nfilter = 50),
                        alphaSet = seq(0.8, 1, 0.1),
                        min_1se = 0, cv.cores = 8)
summary(etaMod)  # AUC 0.7647, BA 0.7313

# ETA model with clinical 14/05/24
data3 <- cbind(data2, clindata3)

set.seed(12, "L'Ecuyer-CMRG")
etaMod <- nestcv.glmnet(y[eta], data3[eta, ],
                        family="binomial",
                        filterFUN = ttest_filter,
                        filter_options = list(nfilter = 50),
                        alphaSet = seq(0.8, 1, 0.1),
                        min_1se = 0, cv.cores = 8)
summary(etaMod)


# TOC model with clinical 29/4/24
trc <- trainControl(method = "cv",
                    number = 10,
                    classProbs = TRUE,
                    savePredictions = "final",
                    summaryFunction = defaultSummary)

# AUC 0.7480; BA 0.7024
data3 <- cbind(data2, clindata3)

tocMod <- nestcv.train(y[toc], data3[toc, ],
                       method="gbm",
                       metric="Kappa",
                       trControl = trc,
                       filterFUN=ttest_filter,
                       filter_options = list(nfilter = 30),
                       balance = "smote",
                       outer_train_predict = T,
                       cv.cores = 8)
summary(tocMod)
tocMod$final_vars

save(rtxMod, etaMod, tocMod,
     file = "/../STRAP models/STRAP_final_models_wclin.Rdata")


## ROC plots (revised with clinical params)

pdf("/../STRAP models/final_rocs_wclin_140524.pdf",
    width=9.5, height=3.2)
op <- par(pty="s", mfrow=c(1, 3), tcl=-0.4)
r1 <- innercv_roc(etaMod)
plot(r1, col="gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Etanercept model (glmnet)", font.main = 1)
lines(etaMod$roc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(etaMod$roc$auc, r1$auc), digits=3)),
       bty="n")

r2 <- innercv_roc(tocMod)
plot(r2, col="gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Tocilizumab model (GBM)", font.main = 1)
lines(tocMod$roc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(tocMod$roc$auc, r2$auc), digits=3)),
       bty="n")

plot(rtx.inroc, col="gold", las=1, mgp=c(2.3, 0.7, 0),
     main="Rituximab model (xgbLinear)", font.main = 1)
lines(rtxroc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(rtxroc$auc, rtx.inroc$auc), digits=3)),
       bty="n")

par(op)
dev.off()


## var_stability plots

vs1 <- plot_var_stability(etaMod, direction = 1, scheme = c("red", "royalblue")) +
  theme(legend.position = "none")
vs2 <- plot_var_stability(tocMod, direction = 1, scheme = c("red", "royalblue")) +
  theme(legend.position = "none")
vs3 <- plot_var_stability(rtxMod, direction = 1, scheme = c("red", "royalblue"),
                          dir_labels = c("Up in Non-responder", "Up in Responder")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0),
         legend.background = element_rect(fill="white",
                                          colour="white"))
  theme(legend.position = c(0.72, 0.25),
        legend.background = element_rect(fill="white",
                                         colour="white"))
# check plot  
vs3
  
pp <- ggarrange(vs1, vs2, vs3, ncol = 3)

ggsave("/../STRAP models/var_stability_wclin_140524.pdf",
       pp, width = 12, height = 7)

## model comparison plots
# trc to handle imbalance

mods <- list(
  glmnet = list(model="glmnet", args=list(family="binomial",
                                          alphaSet = seq(0.7, 1, 0.1),
                                          min_1se = 0)),
  rf = list(model="rf", args=list(filterFUN = "ttest_filter",
                                  filter_options = list(nfilter = 40),
                                  balance="smote")),
  gbm = list(model="gbm", args=list(metric="Kappa",
                                    trControl = trc,
                                    filterFUN = "ttest_filter",
                                    filter_options = list(nfilter = 40))),
  mda = list(model="mda", args=list(metric="Kappa",
                                    trControl = trc,
                                    filterFUN = "ttest_filter",
                                    filter_options = list(nfilter = 40))),
  svmRadial = list(model="svmRadial", args=list(metric="Kappa",
                                                trControl = trc,
                                                filterFUN = "ttest_filter",
                                                filter_options = list(nfilter = 40),
                                                balance="smote")),
  svmPoly = list(model="svmPoly", args=list(metric="Kappa",
                                            trControl = trc,
                                            filterFUN = "ttest_filter",
                                            filter_options = list(nfilter = 40),
                                            balance="smote")),
  xgbTree = list(model="xgbTree", args=list(metric="Kappa",
                                            trControl = trc,
                                            filterFUN = "ttest_filter",
                                            filter_options = list(nfilter = 40))),
  xgbLinear = list(model="xgbLinear", args=list(metric="Kappa",
                                                trControl = trc,
                                                filterFUN = "ttest_filter",
                                                filter_options = list(nfilter = 40)))
)

# ETA
set.seed(12, "L'Ecuyer-CMRG")
resEta <- compare_mod(y[eta], data2[eta, ], mods,
                      repeats = 25, n_outer_folds = 8, rep.cores = 9)

# TOC
resToc <- compare_mod(y[toc], data2[toc, ], mods,
                      repeats = 25, n_outer_folds = 8, rep.cores = 9)

p1 <- plot_modcompare(resEta) + xlab("") + ggtitle("Etanercept")
p2 <- plot_modcompare(resToc) + xlab("") + ggtitle("Tocilizumab")
p3 <- plot_modcompare(resRtx, "Rsquared") + xlab("") + ggtitle("Rituximab") +
  ylab(expr(R^2))

pp <- ggarrange(p1, p2, p3, ncol = 3)

ggsave("/../STRAP models/final_compare_mods_191123.pdf",
       pp, width = 9, height = 3.5)

# clinical only
mods_clin <- list(
  glmnet = list(model="glmnet", args=list(family="binomial")),
  rf = list(model="rf"),
  gbm = list(model="gbm"),
  mda = list(model="mda"),
  svmRadial = list(model="svmRadial"),
  svmPoly = list(model="svmPoly"),
  xgbTree = list(model="xgbTree"),
  xgbLinear = list(model="xgbLinear")
)

clinset5 <- c("Number.of.Tender.Joints", "Number.of.Swollen.Joints", "ESR", "CRP",
              "CCP", "RF", "Arthritis.Activity", "HAQ.Score", "Age", "Gender")
clindata4 <- metadata[, clinset5]

# check individual model
fit <- nestcv.glmnet(y[toc], clindata4[toc, ], family = "binomial", alphaSet = 1)
cx <- checkxy(y[toc], clindata4[toc, ])

# compare models
clinToc <- compare_mod(y[toc], clindata4[toc, ], mods_clin,
                       repeats = 25, n_outer_folds = 8, rep.cores = 13)

clinEta <- compare_mod(y[eta], clindata2[eta, ], mods_clin,
                       repeats = 25, n_outer_folds = 8, rep.cores = 13)
clinRtx <- compare_mod(y[rtx], clindata4[rtx, ], mods_clin,
                       repeats = 25, n_outer_folds = 8, rep.cores = 13)

max(sapply(clinEta$roc, pROC::auc))
max(sapply(clinToc$roc, pROC::auc))
max(sapply(clinRtx$roc, pROC::auc))

p1 <- plot_modcompare(clinEta) + xlab("") + ggtitle("Etanercept")
p2 <- plot_modcompare(clinToc) + xlab("") + ggtitle("Tocilizumab")
p3 <- plot_modcompare(clinRtx) + xlab("") + ggtitle("Rituximab")

pp <- ggarrange(p1, p2, p3, ncol = 3)

ggsave("/../STRAP models/clinical_compare_mods_290424.pdf",
       pp, width = 9, height = 3.5)

clinEta_rocs <- sapply(clinEta$roc, auc)  # xgbLinear
clinToc_rocs <- sapply(clinToc$roc, auc)  # glmnet
clinRtx_rocs <- sapply(clinRtx$roc, auc)  # xgbTree

clinEta_roc <- clinEta$roc[["xgbLinear"]]
clinToc_roc <- clinToc$roc[["glmnet"]]
clinRtx_roc <- clinRtx$roc[["xgbTree"]]

clinEta_roc <- clinEta$roc[["glmnet"]]
clinToc_roc <- clinToc$roc[["gbm"]]
clinRtx_roc <- clinRtx$roc[["xgbLinear"]]

pdf("/../STRAP models/clin_rocs_071124.pdf",
    width=9.5, height=3.2)
op <- par(pty="s", mfrow=c(1, 3), tcl=-0.4)
plot(clinEta_roc, col="green3", las=1, mgp=c(2.3, 0.7, 0),
     main="Etanercept clinical model (glmnet)", font.main = 1)
legend("bottomright", lwd = 2, col = "green3",
       legend=paste("Repeat nested CV AUC", format(clinEta_roc$auc, digits=3)),
       bty="n")
plot(clinToc_roc, col="green3", las=1, mgp=c(2.3, 0.7, 0),
     main="Tocilizumab clinical model (gbm)", font.main = 1)
legend("bottomright", lwd = 2, col = "green3",
       legend=paste("Repeat nested CV AUC",
                    format(clinToc_roc$auc, digits=3)),
       bty="n")
plot(clinRtx_roc, col="green3", las=1, mgp=c(2.3, 0.7, 0),
     main="Rituximab clinical model (xgbLinear)", font.main = 1)
legend("bottomright", lwd = 2, col = "green3",
       legend=paste("Repeat nested CV AUC",
                    format(clinRtx_roc$auc, digits=3)),
       bty="n")
par(op)
dev.off()

# setup clinical models to test in R4RA
set.seed(12, "L'Ecuyer-CMRG")
tocModClin <- nestcv.glmnet(y[toc], clindata4[toc, ],
                            family = "binomial", cv.cores = 8)
summary(tocModClin)

true_y <- meta_RNA_samples$TargetDAS28.ESR.V7[ind]
true_y <- factor(true_y, labels = c("Non.responder", "Responder"))

ind <- meta_RNA_samples$Cohort == "R4RA" & meta_RNA_samples$Visit == 3 &
  meta_RNA_samples$Treatment == "Tocilizumab"

r4ra_toc_dat2 <- vst_all_comBat_clin[ind, ]

r4ra_toc_preds_CL <- as.vector(predict(tocModClin, r4ra_toc_dat2))
r4ra_toc_roc_CL <- pROC::roc(true_y, r4ra_toc_preds_CL, direction = "<")
pROC::auc(r4ra_toc_roc_CL)

save.image("/../STRAP_more_models_140524.RData")

save(etaMod, tocMod, rtxMod, file = "/../STRAP_final_models_140524.RData")

# check seropositivity table
t1 <- table(metadata$ACR.response.status.V7[eta], metadata$RF_status[eta])
R <- t1["Responder", ]
paste0(R, " (", round(R / colSums(t1) *100), "%)")
format(chisq.test(t1)$p.value, digits = 2)

