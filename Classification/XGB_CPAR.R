library(arules)
library(arulesCBA)
library(caTools)
library(MASS)
library(data.table)
library(tidyverse)
library(pROC)
library(PRROC)

source("functions_CPAR_classification.R")

D <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)
                     [which(present_proportion >= 0.1),
                       sampleMetadata_disease
                       [which(sampleMetadata_disease$study_condition == "IBD"),]$sample_id]))

D_false <- -1*D
colnames(D_false) <- paste("-",colnames(D_false))
D <- cbind(D,D_false)
ntaxa <- ncol(D)/2
D$class <- "D"

H <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata$sample_id]))

H_false <- -1*H
colnames(H_false) <- paste("-",colnames(H_false))
H <- cbind(H,H_false)
H$class <- "H"

test_classification <- rbind(H,D)
class <- as.factor(test_classification$class)

set.seed(11)
sample_CPAR_80 <- sample(1:nrow(test_classification),nrow(test_classification)*0.5)
test_20 <- test_classification[-sample_CPAR_80,]
train_80 <- test_classification[sample_CPAR_80,]


train_80_CPAR <- as.data.frame(apply(train_80[,1:ntaxa],
                                     2,function(i) as.logical(as.integer(unlist(i)))))

train_80_CPAR$class <- class[sample_CPAR_80]


set.seed(11)
library(caret)
Folds_CPAR <- createFolds(1:nrow(train_80_CPAR), k = 10)

for (i in 1:10){
  train <- train_80_CPAR[unlist(Folds_CPAR[-i]),]
  #test <- CPAR_20[unlist(Folds_CPAR),]
  cl <- CPAR(class ~ ., train)
  if (i == 1){
    ruleset <- cl[["rules"]]
  }
  else{
    ruleset <- c(ruleset,cl[["rules"]])
  }
  
}

itemset <- CPAR_feature_RA(ruleset)



library(reticulate)



set.seed(1111)

library(caret)
Folds <- createFolds(1:nrow(test_20), k = 10)
accuracy_XGB <- c()
F1_XGB <- c()
auc_XGB <- c()
auprc_XGB <- c()

for (i in 1:10){
  
  train <- rbind(test_20[unlist(Folds[-i]),itemset],train_80[,itemset])
  test <- test_20[unlist(Folds[i]),itemset]
  

  class_train <- as.factor(c(test_20[unlist(Folds[-i]),ncol(test_20)],train_80$class))
  
  XGB <- reticulate :: import("xgboost", convert = FALSE)
  
  xgb_model = XGB$XGBClassifier(objective="binary:logistic")
  classifier_XGB = xgb_model$fit(as.matrix(train), as.matrix(as.numeric(class_train)-1))
  
  y_pred = as.vector(xgb_model$predict(as.matrix(test)))
  y_pred_prob = as.matrix(xgb_model$predict_proba(as.matrix(test)))
  
  
  eval_XGB <- confusionMatrix(as.factor(y_pred), 
                              as.factor(as.numeric(as.factor(as.factor(test_20[unlist(Folds[i]),"class"])))-1), 
                              mode = "everything", positive="0")
  accuracy_XGB <- c(accuracy_XGB,eval_XGB$overall["Accuracy"])
  F1_XGB <- c(F1_XGB,eval_XGB$byClass["F1"])
  roc_XGB <- roc(as.numeric(as.factor(test_20[unlist(Folds[i]),"class"]))-1,
                 y_pred_prob[,2],levels = c("0","1"),direction = "<")
  
  prroc_XGB <- pr.curve(scores.class0 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "D"),1],
                        scores.class1 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "H"),1],curve=T)
  
  auc_XGB <- c(auc_XGB,roc_XGB$auc)
  auprc_XGB <- c(auprc_XGB,prroc_XGB$auc.integral)
  
  
}

mean(accuracy_XGB)
mean(F1_XGB)
mean(auc_XGB)
mean(auprc_XGB)


length(itemset)

IBD_XGB <- rbind(accuracy_XGB,F1_XGB,auc_XGB,auprc_XGB)