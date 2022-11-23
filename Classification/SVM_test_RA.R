# Installing package
#install.packages("caTools")       # For sampling the dataset
#install.packages("randomForest")  # For implementing random forest algorithm

# Loading package
library(caTools)
library(e1071)
library(pROC)
library(PRROC)

library(data.table)
CRC <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata_disease[which(sampleMetadata_disease$study_condition == "T2D"),]$sample_id]))
CRC_false <- -1*CRC
colnames(CRC_false) <- paste("-",colnames(CRC_false))
CRC <- cbind(CRC,CRC_false)
CRC$class <- "CRC"

H <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata$sample_id]))
H_false <- -1*H
colnames(H_false) <- paste("-",colnames(H_false))
H <- cbind(H,H_false)
H$class <- "H"

test_classification <- rbind(CRC,H)
class <- as.factor(test_classification$class)

test_classification$class <- class

set.seed(6666)
library(caret)
Folds <- createFolds(1:nrow(test_classification), k = 10)
accuracy_SVM <- c()
F1_SVM <- c()
auc_SVM <- c()
auprc_SVM <- c()
for (i in 1:10){
  train <- test_classification[unlist(Folds[-i]),]
  test <- test_classification[unlist(Folds[i]),]
  classifier_SVM = svm(formula = class  ~ .,
                       data = train,
                       type = 'C-classification',
                       kernel = 'linear',
                       probability=TRUE)

  y_pred = predict(classifier_SVM, newdata = test[,-ncol(test_classification)])
  y_pred_prob = attr(predict(classifier_SVM, newdata = test[-ncol(test_classification)], probability=TRUE),"probabilities")
  
  eval_SVM <- confusionMatrix(as.factor(as.numeric(y_pred)), as.factor(as.numeric(as.factor(test[,ncol(test_classification)]))), 
                             mode = "everything", positive="1")
  accuracy_SVM <- c(accuracy_SVM,eval_SVM$overall["Accuracy"])
  F1_SVM <- c(F1_SVM,eval_SVM$byClass["F1"])
  roc_SVM <- roc(as.numeric(test[,ncol(test_classification)])-1,
                y_pred_prob[,"H"],levels = c("0","1"),direction = "<")
  
  prroc_SVM <- pr.curve(scores.class0 = y_pred_prob[which(test[,ncol(test_classification)] == "CRC"),"CRC"],
                       scores.class1 = y_pred_prob[which(test[,ncol(test_classification)] == "H"),"CRC"],curve=T)
  
  auc_SVM <- c(auc_SVM,roc_SVM$auc)
  auprc_SVM <- c(auprc_SVM,prroc_SVM$auc.integral)
}

mean(accuracy_SVM)
mean(F1_SVM)
mean(auc_SVM)
mean(auprc_SVM)

performace_list_T2D_SVM_RA <- list("T2D_Acc_SVM" = accuracy_SVM,"T2D_F1_SVM" = F1_SVM,
                                   "T2D_AUC_SVM" = auc_SVM,"T2D_AUPRC_SVM" = auprc_SVM)

fwrite(performace_list_T2D_SVM_RA,"performace_list_T2D_SVM_RA.txt")
