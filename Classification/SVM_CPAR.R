library(e1071)
library(arules)
library(arulesCBA)
library(caTools)
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


set.seed(1111)
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


set.seed(1111)

library(caret)
Folds <- createFolds(1:nrow(test_20), k = 10)
accuracy_SVM <- c()
F1_SVM <- c()
auc_SVM <- c()
auprc_SVM <- c()

for (i in 1:10){
  
  train <- rbind(test_20[unlist(Folds[-i]),itemset],train_80[,itemset])
  test <- test_20[unlist(Folds[i]),itemset]
  
  
  class_train <- as.factor(c(test_20[unlist(Folds[-i]),ncol(test_20)],train_80$class))
  
  classifier_SVM = svm(formula = class_train  ~ .,
                       data = cbind(train,class_train),
                       type = 'C-classification',
                       kernel = 'linear',
                       probability=TRUE)
  
  y_pred = predict(classifier_SVM, newdata = test)
  y_pred_prob = attr(predict(classifier_SVM, newdata = test, probability=TRUE),"probabilities")
  
  
  eval_SVM <- confusionMatrix(as.factor(as.numeric(y_pred)), 
                              as.factor(as.numeric(as.factor(as.factor(test_20[unlist(Folds[i]),"class"])))), 
                              mode = "everything", positive="1")
  accuracy_SVM <- c(accuracy_SVM,eval_SVM$overall["Accuracy"])
  F1_SVM <- c(F1_SVM,eval_SVM$byClass["F1"])
  roc_SVM <- roc(as.numeric(as.factor(test_20[unlist(Folds[i]),"class"]))-1,
                 y_pred_prob[,"H"],levels = c("0","1"),direction = "<")
  
  prroc_SVM <- pr.curve(scores.class0 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "D"),"D"],
                        scores.class1 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "H"),"D"],curve=T)
  
  auc_SVM <- c(auc_SVM,roc_SVM$auc)
  auprc_SVM <- c(auprc_SVM,prroc_SVM$auc.integral)
  
  classifier_SVM_F = svm(formula = class_train  ~ .,
                         data = cbind(train_F,class_train),
                         type = 'C-classification',
                         kernel = 'linear',
                         probability=TRUE)
  
  
}

mean(accuracy_SVM)
mean(F1_SVM)
mean(auc_SVM)
mean(auprc_SVM)


length(itemset)

IBD_SVM <- rbind(accuracy_SVM,F1_SVM,auc_SVM,auprc_SVM)

