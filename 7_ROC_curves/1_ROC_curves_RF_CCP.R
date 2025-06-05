#### ROC curves 
# load packages and data
library(dplyr)
library(purrr)
library(pROC)
library(caret)
library(cowplot)
library(ggplot2)
library(gridExtra)
load('data.RData')
r4ra.meta <- readRDS("./R4RA/R4RAmetadata.rds")
metadata <- readRDS('./STRAP/STRAPmetadata.rds')
metadata$ACPApos <- as.numeric(metadata$ACPApos)
metadata$ACR.response.status <- as.factor(metadata$ACR.response.status)
levels(metadata$ACR.response.status) <- c("Non-Responder", "Responder")
metadata <- metadata %>% filter(Response_criteria != "NA")
metadata <- metadata %>% rename("Response_criteria" = "ACR.response.status")
# prepare R4RA data 
r4ra.meta$TargetDAS28.ESR <- as.factor(r4ra.meta$TargetDAS28.ESR)
levels(r4ra.meta$TargetDAS28.ESR) <- c("Non.Responder", "Responder")
r4ra.meta <- r4ra.meta %>% filter(TargetDAS28.ESR != "NA")
r4ra.meta <- r4ra.meta %>% rename("Response_criteria" = "TargetDAS28.ESR")
treat <- c("Rituximab", "Tocilizumab", "Etanercept")
# function for STRAP ROC curve
roc_models_strap <- function(meta, drug) {
  meta <- meta %>% filter(Randomised.Medication == drug)
  model <- map(c("CCP", "RF"), function(x) {
    rocdrug <- roc(meta$Response_criteria ~ meta[, x], direction = "<")
    predicted_classes <- ifelse(rocdrug$original.predictor > 0.5, "Responder", "Non-Responder")
    cm <- confusionMatrix(table(predicted_classes, meta$Response_criteria))
    write.csv(cm$table, file = paste0("./", drug, "_", x, "_ROC_CFM.csv"))
    # Calculate accuracy and balanced accuracy 
    accuracy <- cm$overall['Accuracy']
    balanced_accuracy <- cm$byClass['Balanced Accuracy']
    
    rocdrug$accuracy <- accuracy
    rocdrug$balanced_accuracy <- balanced_accuracy
    
    return(rocdrug)
  })
  names(model) <- c("CCP", "RF")
  return(model)
}
# function for R4RA ROC curve
roc_models_r4ra <- function(meta, drug) {
  meta <- meta %>% filter(Randomised.Medication == drug)
  model <- map(c("CCP", "RF"), function(x) {
    rocdrug <- roc(meta$Response_criteria ~ meta[, x], direction = "<")
    predicted_classes <- ifelse(rocdrug$original.predictor > 0.5, "Responder", "Non.Responder")
    # Ensure 'predicted_classes' has both 'Responder' and 'Non.Responder'
    predicted_classes <- factor(predicted_classes, levels = c("Responder", "Non.Responder"))
    meta$Response_criteria <- factor(meta$Response_criteria, levels = c("Responder", "Non.Responder"))
    
    # Create a table with the required levels
    conf_matrix_table <- table(predicted_classes, meta$Response_criteria)
    
    # Check if any level is missing and add it with count zero
    if (nrow(conf_matrix_table) != 2 || ncol(conf_matrix_table) != 2) {
      missing_levels <- setdiff(c("Responder", "Non.Responder"), rownames(conf_matrix_table))
      for (lvl in missing_levels) {
        conf_matrix_table <- rbind(conf_matrix_table, setNames(rep(0, ncol(conf_matrix_table)), colnames(conf_matrix_table)))
        rownames(conf_matrix_table)[nrow(conf_matrix_table)] <- lvl
      }
      missing_levels <- setdiff(c("Responder", "Non.Responder"), colnames(conf_matrix_table))
      for (lvl in missing_levels) {
        conf_matrix_table <- cbind(conf_matrix_table, setNames(rep(0, nrow(conf_matrix_table)), rownames(conf_matrix_table)))
        colnames(conf_matrix_table)[ncol(conf_matrix_table)] <- lvl
      }
    }
    
    # Generate the confusion matrix
    cm <- confusionMatrix(conf_matrix_table)
    
    write.csv(cm$table, file = paste0("./", drug, "_", x, "_ROC_CFM_R4RA.csv"))
    # Calculate accuracy and balanced accuracy 
    accuracy <- cm$overall['Accuracy']
    balanced_accuracy <- cm$byClass['Balanced Accuracy']
    
    rocdrug$accuracy <- accuracy
    rocdrug$balanced_accuracy <- balanced_accuracy
    
    return(rocdrug)
  })
  names(model) <- c("CCP", "RF")
  return(model)
}
# create models
models_strap <- map(treat, ~roc_models_strap(meta = metadata,  .x))
names(models_strap) <- treat
models_r4ra <- map(treat[1:2], ~roc_models_r4ra(meta = r4ra.meta,  .x))
names(models_r4ra) <- treat[1:2]

## plot

# Define the plotting function
plot_roc_base <- function(roc_model, title, additional_roc = NULL) {
  plot.roc(roc_model, main = title, col="green3", las=1, mgp=c(2.3, 0.7, 0), font.main = 1,
           cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.4)  # Adjust the cex_value directly
  
  # Plot additional ROC curve for RTX and ETA
  if (!is.null(additional_roc)) {
    lines.roc(additional_roc, col = "purple", lty = 2)  # Add purple dashed line for additional ROC curve
  }
  
  legend("bottomright", 
         legend = c(
           paste("STRAP AUC: ", sprintf("%.2f", roc_model$auc))
         ), 
         col = "green3", bty = "n", cex = 1.4)
  
  if (!is.null(additional_roc)) {
    legend("bottomright", 
           legend = c(
             paste("R4RA AUC: ", sprintf("%.2f", additional_roc$auc))
           ), 
           col = c("green3", "purple"), bty = "n", cex = 1.4, inset=c(0,0.05))
  }
}

# Plot each ROC curve
pdf("roccurves_CCPandRF_R4RA_STRAP.pdf", width = 14, height = 10, paper = "special")
par(mfrow = c(2, 3))

plot_roc_base(models_strap$Etanercept$CCP, "etanercept anti-CCP titre")
plot_roc_base(models_strap$Tocilizumab$CCP, "tocilizumab anti-CCP titre", models_r4ra$Tocilizumab$CCP)
plot_roc_base(models_strap$Rituximab$CCP, "rituximab anti-CCP titre", models_r4ra$Rituximab$CCP)
plot_roc_base(models_strap$Etanercept$RF, "etanercept RF titre")
plot_roc_base(models_strap$Tocilizumab$RF, "tocilizumab RF titre", models_r4ra$Tocilizumab$RF)
plot_roc_base(models_strap$Rituximab$RF, "rituximab RF titre", models_r4ra$Rituximab$RF)

# Close the PDF device
dev.off()