# Installing package
#install.packages("caTools")       # For sampling the dataset
#install.packages("randomForest")  # For implementing random forest algorithm

# Loading package
library(caTools)
library(randomForest)
library(pROC)
library(PRROC)

library(data.table)
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
accuracy_RF <- c()
F1_RF <- c()
auc_RF <- c()
auprc_RF <- c()
for (i in 1:10){
  train <- test_classification[unlist(Folds[-i]),]
  test <- test_classification[unlist(Folds[i]),]
  classifier_RF = randomForest(x = train[,-ncol(test_classification)],
                               y = as.factor(train$class),
                               ntree = 500)
  y_pred = predict(classifier_RF, newdata = test[,-ncol(test_classification)])
  y_pred_prob = predict(classifier_RF, newdata = test[-ncol(test_classification)], type = "prob")
  
  eval_RF <- confusionMatrix(as.factor(as.numeric(y_pred)), as.factor(as.numeric(as.factor(test[,ncol(test_classification)]))), 
                             mode = "everything", positive="1")
  accuracy_RF <- c(accuracy_RF,eval_RF$overall["Accuracy"])
  F1_RF <- c(F1_RF,eval_RF$byClass["F1"])
  roc_RF <- roc(as.numeric(test[,ncol(test_classification)])-1,
                y_pred_prob[,"H"],levels = c("0","1"),direction = "<")
  #roc_RF <- roc(as.numeric(as.factor(test[,ncol(test_classification)])),
                #as.numeric(y_pred),levels = c("2","1"),direction = ">")
  prroc_RF <- pr.curve(scores.class0 = y_pred_prob[which(test[,ncol(test_classification)] == "CRC"),"CRC"],
                       scores.class1 = y_pred_prob[which(test[,ncol(test_classification)] == "H"),"CRC"],curve=T)
  #prroc_RF <- pr.curve(scores.class0 = ifelse(y_pred[which(test[,ncol(test_classification)] == "CRC")] == "CRC",1,0),
                       #scores.class1 = ifelse(y_pred[which(test[,ncol(test_classification)] == "H")] == "CRC",1,0),curve=T)
  auc_RF <- c(auc_RF,roc_RF$auc)
  auprc_RF <- c(auprc_RF,prroc_RF$auc.integral)
}

mean(accuracy_RF)
mean(F1_RF)
mean(auc_RF)
mean(auprc_RF)

performace_list_CRC_RF_RA <- list("CRC_Acc_RF_20" = accuracy_RF,"CRC_F1_RF_20" = F1_RF,"CRC_AUC_RF_20" = auc_RF,"CRC_AUPRC_RF_20" = auprc_RF)
