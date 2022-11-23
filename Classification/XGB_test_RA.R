

library(reticulate)
library(MASS)
library(data.table)
library(tidyverse)
library(pROC)
library(PRROC)

CRC <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata_disease[which(sampleMetadata_disease$study_condition == "CRC"),]$sample_id]))
#CRC_false <- -1*CRC
#colnames(CRC_false) <- paste("-",colnames(CRC_false))
#CRC <- cbind(CRC,CRC_false)
CRC$class <- "CRC"

H <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata$sample_id]))
#H_false <- -1*H
#colnames(H_false) <- paste("-",colnames(H_false))
#H <- cbind(H,H_false)
H$class <- "H"

test_classification <- rbind(CRC,H)
class <- as.factor(test_classification$class)

test_classification$class <- class

set.seed(6666)
library(caret)
Folds <- createFolds(1:nrow(test_classification), k = 10)
accuracy_XGB <- c()
F1_XGB <- c()
auc_XGB <- c()
auprc_XGB <- c()
for (i in 1:10){
  train <- test_classification[unlist(Folds[-i]),]
  test <- test_classification[unlist(Folds[i]),]
  
  
  XGB <- reticulate :: import("xgboost", convert = FALSE)
  
  xgb_model = XGB$XGBClassifier(objective="binary:logistic")
  classifier_XGB = xgb_model$fit(as.matrix(train[,-ncol(train)]), as.matrix(as.numeric(train[,ncol(train)])-1))
  
  y_pred = as.vector(xgb_model$predict(as.matrix(test[,-ncol(test)])))
  y_pred_prob = as.matrix(xgb_model$predict_proba(as.matrix(test[,-ncol(test)])))
  
  
  eval_XGB <- confusionMatrix(as.factor(y_pred), 
                              as.factor(as.numeric(as.factor(as.factor(test[,"class"])))-1), 
                              mode = "everything", positive="0")
  accuracy_XGB <- c(accuracy_XGB,eval_XGB$overall["Accuracy"])
  F1_XGB <- c(F1_XGB,eval_XGB$byClass["F1"])
  roc_XGB <- roc(as.numeric(as.factor(test[,"class"]))-1,
                 y_pred_prob[,2],levels = c("0","1"),direction = "<")
  
  prroc_XGB <- pr.curve(scores.class0 = y_pred_prob[which(as.factor(test[,"class"]) == "CRC"),1],
                        scores.class1 = y_pred_prob[which(as.factor(test[,"class"]) == "H"),1],curve=T)
  
  auc_XGB <- c(auc_XGB,roc_XGB$auc)
  auprc_XGB <- c(auprc_XGB,prroc_XGB$auc.integral)
}

mean(accuracy_XGB)
mean(F1_XGB)
mean(auc_XGB)
mean(auprc_XGB)

performace_list_CRC_XGB_RA <- list("CRC_Acc_XGB" = accuracy_XGB,"CRC_F1_XGB" = F1_XGB,
                                   "CRC_AUC_XGB" = auc_XGB,"CRC_AUPRC_XGB" = auprc_XGB)

fwrite(performace_list_CRC_XGB_RA,"performace_list_CRC_XGB_RA.txt")

