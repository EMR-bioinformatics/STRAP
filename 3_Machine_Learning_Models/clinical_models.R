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
library(nestedcv) 
library(pROC)
library(SuperLearner)

library(parallel)

RhpcBLASctl::omp_set_num_threads(1L) #switch off multithreading for xgboost to work 
source("compare_mod.R")

meta_BL <- readRDS("R4RA_STRAP_meta.rds")

meta_BL$DAS.ESR.status.V7 <- factor(meta_BL$DAS.ESR.status.V7,
                                    levels = c("High.DA",
                                               "Moderate.DA",
                                               "Low.DA",
                                               "Remission"))

####Rituximab####

STRAP_BL <- meta_BL[grepl("STRAP", meta_BL$SeqID), ]

rtx <- STRAP_BL$Treatment == "Rituximab"
toc <- STRAP_BL$Treatment == "Tocilizumab"

# original
oy <- as.vector(STRAP_BL$TargetDAS28.ESR.V7)
oy <- factor(make.names(oy))
names(oy) <- STRAP_BL$SeqID

# use 4 level DAS28.ESR.status

y <- STRAP_BL$DAS.ESR.status.V7
y <- (as.numeric(y) - 1) / 3

clinset2 <- c("sqTJ", "sqSJ", "sqESR", "lnCRP", "Arthritis.Activity", "CCP", "RF")
#rename for ease/later plotting
colnames(STRAP_BL)[colnames(STRAP_BL) == "logCRP"] <- "lnCRP"
rownames(STRAP_BL) <- STRAP_BL$SeqID
clindata2 <- STRAP_BL[ , clinset2]

xgbl <- list(
  xgbLinear = list(model="xgbLinear")
)

rtx4 <- compare_mod(y[rtx], clindata2[rtx, ], xgbl,
                    repeats = 25, n_outer_folds = 10, rep.cores = 8, 
                    return_fits = T)

#saveRDS(rtx4, "Rit_4levels_xgbl.rds")

#combined auc for outerfolds

rtx4_outer <- data.frame()

for(i in 1:25){
  for(j in 1:10){
    temp <- rtx4$fits$xgbLinear[[i]]$outer_result[[j]]$preds
    temp$oy <- oy[match(rownames(temp), names(oy))]
    temp$rep <- paste(i, "-", j)
    rtx4_outer <- rbind(rtx4_outer, temp)
  }
}

roc(rtx4_outer$oy, rtx4_outer$predy)
#AUC 0.642 for combined train

rtx4_outer$rep_num <- as.numeric(factor(rtx4_outer$rep))

#test in R4RA
R4RA_BL <- meta_BL[grepl("R4RA", meta_BL$SeqID), ]

rtx_test <- R4RA_BL$Treatment == "Rituximab"
toc_test <- R4RA_BL$Treatment == "Tocilizumab"

colnames(R4RA_BL)[colnames(R4RA_BL) == "logCRP"] <- "lnCRP"
rownames(R4RA_BL) <- R4RA_BL$SeqID
clindata_test <- R4RA_BL[ , clinset2]

# original
true_y <- as.vector(R4RA_BL$TargetDAS28.ESR.V7)
true_y <- factor(make.names(true_y))

rtx4_test <- as.numeric()

for(i in 1:25){
  for(j in 1:10){
    train <- rtx4$fits$xgbLinear[[i]]$outer_result[[j]]$fit
    test.df <- data.frame("Pred" = predict(train, clindata_test[rtx_test, ]),
                          "true_y" = true_y[rtx_test],
                          "rep" = paste(i, "-", j))
    rtx4_test <- rbind(rtx4_test, test.df)
  }
}  

roc(rtx4_test$true_y, rtx4_test$Pred)
#AUC 0.6532 for combined test

####Tocilizumab####

gbm <- list(
  gbm = list(model="gbm")
)

tocMod <- compare_mod(oy[toc], clindata2[toc, ], gbm,
                      repeats = 25, n_outer_folds = 10, rep.cores = 8, 
                      return_fits = T)

#saveRDS(tocMod, "Toc_gbm.rds")

toc_outer <- data.frame()

for(i in 1:25){
  for(j in 1:10){
    temp <- tocMod$fits$gbm[[i]]$outer_result[[j]]$preds
    temp$rep <- paste(i, "-", j)
    toc_outer <- rbind(toc_outer, temp)
  }
}

roc(toc_outer$testy, toc_outer$predyp)
#AUC 0.6732 for combined train

#test in R4RA
toc_testres <- as.numeric()

for(i in 1:25){
  for(j in 1:10){
    train <- tocMod$fits$gbm[[i]]$outer_result[[j]]$fit
    test.df <- data.frame("Pred" = predict(train, clindata_test[toc_test, ], type = "prob")[ ,2],
                          "true_y" = true_y[toc_test],
                          "rep" = paste(i, "-", j))
    toc_testres <- rbind(toc_testres, test.df)
  }
}  

roc(toc_testres$true_y, toc_testres$Pred)
#AUC 0.6242 for combined test

####exporting AUC plots####

pdf("clinical_rocs.pdf",
    width=9.5, height=3.2)
op <- par(pty="s", mfrow=c(1, 3), tcl=-0.4)

plot.new()

plot(roc(toc_outer$testy, toc_outer$predyp), col = "green3",
     las = 1, mgp=c(2.3, 0.7, 0),
     main="Tocilizumab model (GBM)", font.main = 1)
lines(roc(toc_testres$true_y, toc_testres$Pred), col = "purple")
legend("bottomright", 
       col = c("green3", "purple"),
       text.col = "black",
       lty = 1,
       legend=paste(c("STRAP nested CV AUC", "Retest in R4RA AUC"),
                    format(c(roc(toc_outer$testy, toc_outer$predyp)$auc,
                             roc(toc_testres$true_y, toc_testres$Pred)$auc),
                           digits=3)),
       bty="n")

plot(roc(rtx4_outer$oy, rtx4_outer$predy), col = "green3",
     las = 1, mgp=c(2.3, 0.7, 0),
     main="Rituximab model (xgbLinear)", font.main = 1)
lines(roc(rtx4_test$true_y, rtx4_test$Pred), col = "purple")
legend("bottomright", 
       col = c("green3", "purple"),
       text.col = "black",
       lty = 1,
       legend=paste(c("STRAP nested CV AUC", "Retest in R4RA AUC"),
                    format(c(roc(rtx4_outer$oy, rtx4_outer$predy)$auc,
                             roc(rtx4_test$true_y, rtx4_test$Pred)$auc),
                           digits=3)),
       bty="n")
par(op)
dev.off()
