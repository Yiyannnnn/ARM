library(MASS)
library(glmnet)
library(data.table)
library(tidyverse)
library(pROC)
library(PRROC)

CRC <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata_disease[which(sampleMetadata_disease$study_condition == "IBD"),]$sample_id]))
CRC$class <- "IBD"

H <- as.data.frame(t(as.data.frame(stool_phyloseq@otu_table)[which(present_proportion >= 0.1),sampleMetadata$sample_id]))
H$class <- "H"

test_classification <- rbind(H,CRC)
class <- as.factor(test_classification$class)
#test_classification <- as.matrix(apply(test_classification[,1:200],2,
                                           #function(i) as.integer(unlist(i))))

#test_classification <- cbind(test_classification,class)


set.seed(1)

library(caret)
Folds <- createFolds(1:nrow(test_classification), k = 10)
accuracy_logistic_lasso <- c()
F1_logistic_lasso <- c()
auc_logistic_lasso <- c()
auprc_logiatic_lasso <- c()
for (i in 1:10){
  train <- test_classification[unlist(Folds[-i]),]
  test <- test_classification[unlist(Folds[i]),]
  grid = 10^seq(10, -2, length = 100)
  
  lasso.mod = glmnet(as.matrix(train[,-ncol(test_classification)]), as.numeric(as.factor(train[,ncol(test_classification)])), alpha = 1, lambda = grid)
  
  cv.out = cv.glmnet(as.matrix(train[,-ncol(test_classification)]), as.numeric(as.factor(train[,ncol(test_classification)])), alpha = 1)
  coef.min = coef(cv.out, s = "lambda.min")
  
  logmodel <- glm(as.numeric(as.factor(class))-1 ~.,family = binomial(link = logit) ,data = as.data.frame(train[,-c(which(coef.min == 0)-1)]))
  
  probabilities = logmodel %>% predict(as.data.frame(test[,-ncol(test_classification)]), type = "response")
  y_pred = as.factor(ifelse(probabilities < 0.5, "H", "IBD"))
  #y_pred <- stats :: predict(logmodel,as.data.frame(test[,-201]))
  
  eval_logistic_lasso <- confusionMatrix(as.factor(y_pred), as.factor(test[,ncol(test_classification)]), mode = "everything", positive="IBD")
  accuracy_logistic_lasso <- c(accuracy_logistic_lasso,eval_logistic_lasso$overall["Accuracy"])
  F1_logistic_lasso <- c(F1_logistic_lasso,eval_logistic_lasso$byClass["F1"])
  roc_logistic_lasso <- roc(as.numeric(as.factor(unlist(as.vector(test[,ncol(test_classification)]))))-1,probabilities,levels = c("0","1"),direction = "<")
  #roc_logistic_lasso <- roc(as.numeric(as.factor(unlist(as.vector(test[,ncol(test_classification)])))),as.numeric(as.factor(y_pred)),levels = c("1","2"),direction = "<")
 
  prroc_logistic_lasso <- pr.curve(scores.class0 = probabilities[which(test[,ncol(test_classification)] == "IBD")],
                                   scores.class1 = probabilities[which(test[,ncol(test_classification)] == "H")],curve=T)
  #prroc_logistic_lasso <- pr.curve(scores.class0 = ifelse(y_pred[which(test[,ncol(test_classification)] == "IBD")] == "IBD",1,0),
                                   #scores.class1 = ifelse(y_pred[which(test[,ncol(test_classification)] == "H")] == "IBD",1,0),curve=T)
  auc_logistic_lasso <- c(auc_logistic_lasso,roc_logistic_lasso$auc)
  auprc_logiatic_lasso <- c(auprc_logiatic_lasso,prroc_logistic_lasso$auc.integral)
}

mean(accuracy_logistic_lasso)
mean(F1_logistic_lasso)
mean(auc_logistic_lasso)
mean(auprc_logiatic_lasso)

performace_list_T2D_logit_lasso_RA <- list("T2D_Acc_logitlasso_RA" = accuracy_logistic_lasso,"T2D_F1_logitlasso_RA" = F1_logistic_lasso,
                                        "T2D_AUC_logitlasso_RA" = auc_logistic_lasso,"T2D_AUPRC_logitlasso_RA" = auprc_logiatic_lasso)


