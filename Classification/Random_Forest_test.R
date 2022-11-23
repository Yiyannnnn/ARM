# Installing package
#install.packages("caTools")       # For sampling the dataset
#install.packages("randomForest")  # For implementing random forest algorithm

# Loading package
library(caTools)
library(randomForest)
library(pROC)
library(PRROC)

library(data.table)
#CRC <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata_disease[which(sampleMetadata_disease$study_condition == "CRC"),]$sample_id]))
CRC <- fread("CPAR_itemset_CRC_20.csv")
#CRC <- as.data.frame(CRC_H_rule_both)
CRC$class <- "CRC"
#H <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata$sample_id]))
H <- fread("CPAR_itemset_H(CRC)_20.csv")
#H <- as.data.frame(H_CRC_rule_both)
H$class <- "H"

test_classification <- rbind(CRC,H)
class <- as.factor(test_classification$class)
#test_classification <- as.data.frame(apply(test_classification[,1:(ncol(test_classification)-1)],2,
                                           #function(i) as.logical(as.integer(unlist(i)))))

test_classification$class <- class
colnames(test_classification) <- c(1:(ncol(test_classification)-1),"class")

test_classification <- test_classification[-sample_CPAR_20,]

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
  classifier_RF = randomForest(x = train[,-c("class")],
                               y = as.factor(train$class),
                               ntree = 500)
  y_pred = predict(classifier_RF, newdata = test[,-c("class")])
  y_pred_prob = predict(classifier_RF, newdata = test[,-c("class")], type = "prob")
  
  eval_RF <- confusionMatrix(as.factor(as.numeric(y_pred)), as.factor(as.numeric(as.factor(test$class))), 
                             mode = "everything", positive="1")
  accuracy_RF <- c(accuracy_RF,eval_RF$overall["Accuracy"])
  F1_RF <- c(F1_RF,eval_RF$byClass["F1"])
  roc_RF <- roc(as.numeric(test$class)-1,
                y_pred_prob[,"H"],levels = c("0","1"),direction = "<")
  #roc_RF <- roc(as.numeric(as.factor(test[,ncol(test_classification)])),
                #as.numeric(y_pred),levels = c("2","1"),direction = ">")
  prroc_RF <- pr.curve(scores.class0 = y_pred_prob[which(test[,c("class")] == "CRC"),"CRC"],
                       scores.class1 = y_pred_prob[which(test[,c("class")] == "H"),"CRC"],curve=T)
  #prroc_RF <- pr.curve(scores.class0 = ifelse(y_pred[which(test[,ncol(test_classification)] == "CRC")] == "CRC",1,0),
                       #scores.class1 = ifelse(y_pred[which(test[,ncol(test_classification)] == "H")] == "CRC",1,0),curve=T)
  auc_RF <- c(auc_RF,roc_RF$auc)
  auprc_RF <- c(auprc_RF,prroc_RF$auc.integral)
}

mean(accuracy_RF)
mean(F1_RF)
mean(auc_RF)
mean(auprc_RF)

performace_list_CRC_RF_20 <- list("CRC_Acc_RF_20" = accuracy_RF,"CRC_F1_RF_20" = F1_RF,"CRC_AUC_RF_20" = auc_RF,"CRC_AUPRC_RF_20" = auprc_RF)

fwrite(performace_list_CRC_RF_20,"performace_list_CRC_RF_20.txt")

accuracy_RF_b <- c()
F1_RF_b <- c()
auc_RF_b <- c()

for (j in 1:10){
  test_classification_b <- test_classification[c(which(test_classification$class == "CRC"),sample(1:nrow(H),nrow(CRC))),]
  random_split <- sample(1:nrow(test_classification_b),0.5*nrow(test_classification_b))
  train <- test_classification_b[random_split,]
  test <- test_classification_b[-random_split,]
  classifier_RF = randomForest(x = train[,-201],
                               y = train$class,
                               ntree = 100)
  y_pred = predict(classifier_RF, newdata = test[-201])
  
  eval_RF_b <- confusionMatrix(as.factor(as.numeric(y_pred)), as.factor(as.numeric(test[,201])), 
                             mode = "everything", positive="1")
  accuracy_RF_b <- c(accuracy_RF_b,eval_RF_b$overall["Accuracy"])
  F1_RF_b <- c(F1_RF_b,eval_RF_b$byClass["F1"])
  roc_RF_b <- roc(as.numeric(test[,201]),
                as.numeric(y_pred),levels = c("2","1"),direction = ">")
  auc_RF_b <- c(auc_RF_b,roc_RF$auc)
}


mean(accuracy_RF_b)
mean(F1_RF_b)
mean(auc_RF_b)
