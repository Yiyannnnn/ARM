source("functions_CPAR_classification.R")

library(MASS)
library(glmnet)
library(data.table)
library(tidyverse)
library(pROC)
library(PRROC)

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


set.seed(111)
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


set.seed(111)

library(caret)
Folds <- createFolds(1:nrow(test_20), k = 10)
accuracy_logistic_lasso <- c()
F1_logistic_lasso <- c()
auc_logistic_lasso <- c()
auprc_logiatic_lasso <- c()


for (i in 1:10){
  train <- rbind(test_20[unlist(Folds[-i]),itemset],train_80[,itemset])
  test <- test_20[unlist(Folds[i]),itemset]
  
  grid = 10^seq(10, -2, length = 100)
  
  
  
  class_train <- c(test_20[unlist(Folds[-i]),ncol(test_20)],train_80$class)
  lasso.mod = glmnet(as.matrix(train), as.numeric(as.factor(class_train)), alpha = 1, lambda = grid)
  
  cv.out = cv.glmnet(as.matrix(train), as.numeric(as.factor(class_train)), alpha = 1)
  coef.min = coef(cv.out, s = "lambda.min")
  
  logmodel <- glm(as.numeric(as.factor(class_train))-1 ~.,family = binomial(link = logit) ,
                  data = as.data.frame(train[,-c(which(coef.min == 0)-1)]))
  
  probabilities = logmodel %>% predict(as.data.frame(test), type = "response")
  y_pred = as.factor(ifelse(probabilities < 0.5, "D", "H"))
  #y_pred <- stats :: predict(logmodel,as.data.frame(test[,-201]))
  
  
  eval_logistic_lasso <- confusionMatrix(as.factor(y_pred), as.factor(test_20[unlist(Folds[i]),"class"]), mode = "everything", positive="D")
  accuracy_logistic_lasso <- c(accuracy_logistic_lasso,eval_logistic_lasso$overall["Accuracy"])
  F1_logistic_lasso <- c(F1_logistic_lasso,eval_logistic_lasso$byClass["F1"])
  roc_logistic_lasso <- roc(as.numeric(as.factor(unlist(as.vector(test_20[unlist(Folds[i]),"class"]))))-1,probabilities,levels = c("0","1"),direction = "<")
  
  prroc_logistic_lasso <- pr.curve(scores.class0 = 1-probabilities[which(test_20[unlist(Folds[i]),"class"] == "D")],
                                   scores.class1 = 1-probabilities[which(test_20[unlist(Folds[i]),"class"] == "H")],curve=T)
  
  auc_logistic_lasso <- c(auc_logistic_lasso,roc_logistic_lasso$auc)
  auprc_logiatic_lasso <- c(auprc_logiatic_lasso,prroc_logistic_lasso$auc.integral)
}

mean(accuracy_logistic_lasso)
mean(F1_logistic_lasso)
mean(auc_logistic_lasso)
mean(auprc_logiatic_lasso)


length(itemset)

IBD_logit <- rbind(accuracy_logistic_lasso,F1_logistic_lasso,auc_logistic_lasso,auprc_logiatic_lasso)

