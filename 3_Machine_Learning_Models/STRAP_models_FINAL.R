library(nestedcv)
library(caret)
library(gbm)
library(ggpubr)
source("/../compare_mod.R")

load("/../STRAP_load_data_230521.RData")
nanostring <- readRDS("/../STRAP/geneAnnot_nCounter_custompanel2.rds")
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
# colnames(data2)[nearZeroVar(data2)]  # check which genes are removed
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
clinset2 <- c("sqTJ", "sqSJ", "sqESR", "lnCRP", "Arthritis.Activity", "CCP", "RF")
clindata2 <- metadata[, clinset2]
correl_filter(y, clindata2, type="full")

nz <- nearZeroVar(t(vstdata))
data3 <- t(vstdata[-nz, ])
cm <- matrixStats::colMedians(data3)
# hist(cm)
table(cm > 8)
data3 <- data3[, cm > 8]  # remove low expressed genes
data3 <- cbind(data3, clindata2)
data3 <- as.matrix(data3)

ind <- removeAC(data3)  # remove pseudogenes
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
set.seed(12, "L'Ecuyer-CMRG")
rtxMod <- nestcv.train(y[rtx], data3[rtx, ],
                       method = "rf",
                       filterFUN = glmnet_filter,
                       filter_options = list(nfilter = 28),
                       cv.cores = 10)
summary(rtxMod)  # RMSE 0.3026  Rsq 0.1779

rtxMod$output$oy <- oy[match(rownames(rtxMod$output), names(oy))]
rtxroc <- roc(rtxMod$output$oy, rtxMod$output$predy, direction="<")
incv <- innercv_preds(rtxMod)
incv$oy <- factor(incv$testy > 0.5, labels = c("Non.responder", "Responder"))
rtx.inroc <- roc(incv$oy, incv$predy, direction = "<")
plot(rtxroc, col = "royalblue", las = 1)
lines(rtx.inroc, col = "gold")
rtxroc$auc  # 0.7208
rtx.inroc$auc  # 0.8525

# rtxMod$output$predclass <- factor(rtxMod$output$predy > 0.5, labels = levels(oy))
# confusionMatrix(rtxMod$output$predclass, rtxMod$output$oy, positive = "Responder")

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
  #                                   filter_options = list(nfilter = 30))),
  svmRadial = list(model="svmRadial", args=list(filterFUN = "glmnet_filter",
                                                filter_options = list(nfilter = 30))),
  svmPoly = list(model="svmPoly", args=list(filterFUN = "glmnet_filter",
                                            filter_options = list(nfilter = 30))),
  xgbTree = list(model="xgbTree", args=list(filterFUN = "glmnet_filter",
                                            filter_options = list(nfilter = 30))),
  xgbLinear = list(model="xgbLinear", args=list(filterFUN = "glmnet_filter",
                                                filter_options = list(nfilter = 30)))
)
# mda sometimes fails (error)

resRtx <- compare_mod(y[rtx], data3[rtx, ], linmods,
                      repeats = 25, n_outer_folds = 8, cv.cores = 8)
plot_modcompare(resRtx, "RMSE")
plot_modcompare(resRtx, "Rsquared")

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

# TOC
trc <- trainControl(method = "cv",
                    number = 10,
                    classProbs = TRUE,
                    savePredictions = "final",
                    summaryFunction = defaultSummary)

tocMod <- nestcv.train(y[toc], data2[toc, ],
                       method="gbm",
                       metric="Kappa",
                       trControl = trc,
                       filterFUN=ttest_filter,
                       filter_options = list(nfilter=35, p_cutoff=0.1),
                       balance = "smote",
                       outer_train_predict = T,
                       cv.cores = 8)
summary(tocMod)

save(rtxMod, etaMod, tocMod,
     file = "/../STRAP models/STRAP_final_models.Rdata")


## ROC plots

pdf("/../STRAP models/final_rocs.pdf",
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
     main="Rituximab model (RF)", font.main = 1)
lines(rtxroc, col="blue")
legend("bottomright", legend=paste(c("Nested CV AUC", "Inner CV AUC"),
                                   format(c(rtxroc$auc, rtx.inroc$auc), digits=3)),
       bty="n")

par(op)
dev.off()


## var_stability plots

vs1 <- plot_var_stability(etaMod, direction = 1) +
  theme(legend.position = "none")
vs2 <- plot_var_stability(tocMod, direction = 1) +
  theme(legend.position = "none")
vs3 <- plot_var_stability(rtxMod, direction = 1,
                          dir_labels = c("Up in Non-responder", "Up in Responder")) +
  # theme(legend.justification=c(1,0), legend.position=c(1,0),
  #       legend.background = element_rect(fill="white",
  #                                        colour="white"))
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill="white",
                                         colour="white"))
pp <- ggarrange(vs1, vs2, vs3, ncol = 3)

ggsave("/users/myles/dropbox/R scripts/STRAP/STRAP models/final var_stability.pdf",
       pp, width = 12, height = 7)

## model comparison plots

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
set.seed(28, "L'Ecuyer-CMRG")
resEta <- compare_mod(y[eta], data2[eta, ], mods,
                      repeats = 25, n_outer_folds = 8, cv.cores = 8)

# TOC
data2 <- cbind(data2, clindata2)

resToc <- compare_mod(y[toc], data2[toc, ], mods,
                      repeats = 25, n_outer_folds = 8, cv.cores = 8)

p1 <- plot_modcompare(resEta) + xlab("") + ggtitle("Etanercept")
p2 <- plot_modcompare(resToc) + xlab("") + ggtitle("Tocilizumab")
p3 <- plot_modcompare(resRtx, "Rsquared") + xlab("") + ggtitle("Rituximab")

pp <- ggarrange(p1, p2, p3, ncol = 3)

ggsave("/../STRAP models/final_compare_mods_v2.pdf",
       pp, width = 9, height = 3.5)



save.image("/../STRAP_final_models_250223.RData")
