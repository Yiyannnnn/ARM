library(data.table)
library(stringr)
rule_H <- fread("univ_rule_H_supp0.7_uncleaned.csv")


lset_H = lapply(rule_H$V2,function(i) unlist(strsplit(str_sub(i,3,-3), ", ")))
lset_H_dt <- suppressWarnings(do.call(rbind.data.frame, lset_H))

rset_H = lapply(rule_H$V3,function(i) unlist(strsplit(str_sub(i,3,-3), ", ")))
rset_H_dt <- suppressWarnings(do.call(rbind.data.frame, rset_H))

set_H <- as.data.table(cbind(lset_H_dt,rset_H_dt))

set_list_H <- list()
set_list_H <- apply(set_H,1,function(i) c(set_list_H,i))
set_list_H <- lapply(set_list_H,function(i) unique(unlist(i)))

freq_itemset_H <- unique(lapply(set_list_H,sort))

names(sort(table(unlist(freq_itemset_H)),decreasing = TRUE))

freq_itemset_H_filter <- lapply(freq_itemset_H,function(i) ifelse(all(grepl("-", unlist(i))),NA,i[which(!grepl("-", unlist(i)))]))

freq_itemset_H_filter <- freq_itemset_H[!is.na(freq_itemset_H_filter)]

library(stringi)
freq_itemset_H_dt <- stri_list2matrix(freq_itemset_H_filter, byrow=TRUE)

library(reshape2)
present_itemset_H <- table(reshape2 :: melt(freq_itemset_H_dt, na.rm=TRUE)[,-2])
present_itemset_H_df <- as.data.frame.matrix(present_itemset_H)

library(RInSp)
NODF_H <- NODF(import.RInSp(present_itemset_H_df))

library(vegan)
NTC_H <- nestedtemp(present_itemset_H_df,niter = 100)
plot(NTC_H,names = c(FALSE,TRUE),kind = c("incidence"))

sort(names(NTC_H$c))

library(ggplot2)
library(RColorBrewer)
library(tidyverse)

pdf("./Nestedness_calculation/H_incidence_matrix_filtered.pdf",width=6, height=10)
as.data.table(present_itemset_H)|> 
  mutate(Var1 = factor(Var1,levels = rev(names(NTC_H$r)))) |> 
  mutate(value = factor(value,levels = names(NTC_H$c))) |>
  mutate(class = case_when(grepl("-", value) ~ paste("abs",N), TRUE ~ paste("pre",N))) |>
  #ggplot(aes(value,Var1, fill = class)) +
  ggplot(aes(value,Var1, fill = N)) +
  scale_x_discrete() +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradientn(colors = c("ivory","#008500")) +
  #scale_fill_manual(values = c("blue","green","blue","yellow")) +
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        legend.position="bottom", 
        text = element_text(size = 15),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0, hjust=-0.05))+
  labs(title = "", x = "species", y = "Frequent Itemset")
dev.off()

OTU_H <- fread("transaction_H_filter_10.csv")
rulelist_H <- rules(
  lhs = lapply(lapply(rule_H$V2,function(i) 
    unlist(strsplit(str_sub(i,3,-3), ", "))),function(j) str_sub(j,2,-2)), 
  
  rhs = lapply(lapply(rule_H$V3,function(i)
    unlist(strsplit(str_sub(i,3,-3), ", "))),function(j) str_sub(j,2,-2)), 
  itemLabels = colnames(OTU_H),
  quality = rule_H[,c(4:6)]
)

inspect(head(rulelist_H, n = 10, by = "lift"))


