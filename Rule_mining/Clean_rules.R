#Clean rules

library(tidyverse)
library(data.table)
Rules_py_H <- fread("python_rules_H_0.1pre_0.7supp.csv")
Rules_py_H <- na.omit(Rules_py_H)

H_lhs = lapply(Rules_py_H$antecedents,function(i) strsplit(i, "frozenset")[[1]][2])
H_rhs = lapply(Rules_py_H$consequents,function(i) strsplit(i, "frozenset")[[1]][2])


Rules_py_H <- cbind(Rules_py_H,unlist(H_lhs),unlist(H_rhs))

Rules_py_H$rule <- paste(Rules_py_H$V2, Rules_py_H$V3)

univ_rule_py_H <- NULL
univ_rule_names_H <- names(which(table(Rules_py_H$rule) >= 4))
for (i in univ_rule_names_H){
  rules <- Rules_py_H[which(Rules_py_H$rule == i),]
  rulei <- as.data.frame(apply(rules[,c(6,7,8)],2,mean))
  rulei <- cbind(rules[1,c(11:13)],t(rulei))
  univ_rule_py_H <- rbind(univ_rule_py_H,rulei)
}

dim(univ_rule_py_H)

#fwrite(univ_rule_py_H,"univ_rule_H_supp0.7_uncleaned.csv")

