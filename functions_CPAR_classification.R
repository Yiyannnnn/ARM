library(arules)
CPAR_feature_RA <- function(ruleset){
  ruleset_df <- DATAFRAME(ruleset)
  ruleset_df_univ <- ruleset_df[which(ruleset_df$LHS %in% names(which(table(ruleset_df$LHS) >=2))),c("LHS","RHS","laplace")]
  #ruleset_df_univ <- ruleset_df_univ[which(ruleset_df_univ$laplace >= 0.9),]
  ruleset_df_univ_H <- unique(ruleset_df_univ[which(ruleset_df_univ$RHS == "{class=H}"),"LHS"])
  ruleset_df_univ_H_laplace <- unlist(lapply(ruleset_df_univ_H,function(i) mean(ruleset_df_univ[which(ruleset_df_univ$LHS == i),"laplace"])))
  ruleset_df_univ_D <- unique(ruleset_df_univ[which(ruleset_df_univ$RHS == "{class=D}"),"LHS"])
  ruleset_df_univ_D_laplace <- unlist(lapply(ruleset_df_univ_D,function(i) mean(ruleset_df_univ[which(ruleset_df_univ$LHS == i),"laplace"])))
  
  CPAR_itemset_H <- lapply(ruleset_df_univ_H,function(i) unlist(strsplit(str_sub(i,2,-2), ",")))
  CPAR_itemset_H <- lapply(CPAR_itemset_H,function(i) unlist(lapply(i,function(j) ifelse(as.logical(unlist(strsplit(j, "="))[2]),
                                                                                         unlist(strsplit(j, "="))[1],
                                                                                         paste("-",unlist(strsplit(j, "=")))[1]))))
  
  CPAR_itemset_D <- lapply(ruleset_df_univ_D,function(i) unlist(strsplit(str_sub(i,2,-2), ",")))
  CPAR_itemset_D <- lapply(CPAR_itemset_D,function(i) unlist(lapply(i,function(j) ifelse(as.logical(unlist(strsplit(j, "="))[2]),
                                                                                         unlist(strsplit(j, "="))[1],
                                                                                         paste("-",unlist(strsplit(j, "=")))[1]))))
  
  CPAR_itemset <- append(CPAR_itemset_D,CPAR_itemset_H)
  CPAR_rule_laplace <- c(ruleset_df_univ_D_laplace,ruleset_df_univ_H_laplace)
  
  itemset <- unique(c(unlist(CPAR_itemset_H),unlist(CPAR_itemset_D)))
  
  
  
  return(itemset)
  
}



