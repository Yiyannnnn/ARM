library(arules)
library(arulesCBA)
library(caTools)
library(randomForest)
library(MASS)
library(glmnet)
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

set.seed(1)
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
accuracy_RF <- c()
F1_RF <- c()
auc_RF <- c()
auprc_RF <- c()

accuracy_RF_F <- c()
F1_RF_F <- c()
auc_RF_F <- c()
auprc_RF_F <- c()
for (i in 1:10){
  
  train <- rbind(test_20[unlist(Folds[-i]),itemset],train_80[,itemset])
  test <- test_20[unlist(Folds[i]),itemset]
  
  
  class_train <- c(test_20[unlist(Folds[-i]),ncol(test_20)],train_80$class)
  
  classifier_RF = randomForest(x = train,
                               y = as.factor(class_train),
                               ntree = 500)
  y_pred = predict(classifier_RF, newdata = test)
  y_pred_prob = predict(classifier_RF, newdata = test, type = "prob")
  
  eval_RF <- confusionMatrix(as.factor(as.numeric(y_pred)), 
                             as.factor(as.numeric(as.factor(as.factor(test_20[unlist(Folds[i]),"class"])))), 
                             mode = "everything", positive="1")
  accuracy_RF <- c(accuracy_RF,eval_RF$overall["Accuracy"])
  F1_RF <- c(F1_RF,eval_RF$byClass["F1"])
  roc_RF <- roc(as.numeric(as.factor(test_20[unlist(Folds[i]),"class"]))-1,
                y_pred_prob[,"H"],levels = c("0","1"),direction = "<")
  
  prroc_RF <- pr.curve(scores.class0 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "D"),"D"],
                       scores.class1 = y_pred_prob[which(as.factor(test_20[unlist(Folds[i]),"class"]) == "H"),"D"],curve=T)
  
  auc_RF <- c(auc_RF,roc_RF$auc)
  auprc_RF <- c(auprc_RF,prroc_RF$auc.integral)
  
  
}

mean(accuracy_RF)
mean(F1_RF)
mean(auc_RF)
mean(auprc_RF)


length(itemset)

IBD_RF <- rbind(accuracy_RF,F1_RF,auc_RF,auprc_RF)