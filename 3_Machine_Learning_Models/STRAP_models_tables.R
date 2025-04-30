library(nestedcv)
library(caret)
library(gbm)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(SuperLearner)
library(pROC)
source("/../compare_mod.R")

load("/../STRAP_more_models_290424.RData")

# make nice table
library(reshape2)
mod_table <- function(x) {
  s <- summary(x)
  m <- dcast(s, model ~ metric, value.var = "mean")
  sem <- dcast(s, model ~ metric, value.var = "sem")
  m2 <- as.matrix(format(m[, -1], digits = 3))
  sem2 <- as.matrix(format(sem[ , -1], digits = 2))
  p <- paste(m2, sem2, sep = " + ")
  p <- matrix(p, nrow = nrow(m2))
  df <- cbind(m[, 1], p)
  colnames(df) <- colnames(m)
  df
}
m1 <- mod_table(resEta)
m2 <- mod_table(resToc)
m3 <- mod_table(resRtx)

write.csv(m1, "/../STRAP models/Tables/resEta table.csv")
write.csv(m2, "/../STRAP models/Tables/resToc table.csv")
write.csv(m3, "/../STRAP models/Tables/resRtx table.csv")


save.image("/../STRAP/STRAP_more_models_140524.RData")

save(etaMod, tocMod, rtxMod, file = "/../STRAP/STRAP_final_models_140524.RData")
