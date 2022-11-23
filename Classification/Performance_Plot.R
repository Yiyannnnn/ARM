#Performance Plot

performace_list <- readRDS("performace_list_final.RData")

pdf("performace_by_dis_eval_IBD.pdf",height = 6, width = 3)
performace_list |> filter(str_detect(variable,"IBD")) |> 
  mutate(method = as.factor(case_when(str_detect(variable,"logistic w/ RA") ~ "logistic w/ RA",str_detect(variable,"logistic w/ CPAR") ~ "logistic w/ CPAR",
                                      str_detect(variable,"RF w/ RA") ~ "RF w/RA", str_detect(variable,"RF w/ CPAR") ~ "RF w/CPAR",
                                      str_detect(variable,"SVM w/ RA") ~ "SVM w/ RA ", str_detect(variable,"SVM w/ CPAR") ~ "SVM w/ CPAR",
                                      str_detect(variable,"XGB w/ RA") ~ "XGB w/ RA",str_detect(variable,"XGB w/ CPAR") ~ "XGB w/ CPAR"))) |>
  mutate(disease = case_when(str_detect(variable,"CRC") ~ "CRC",str_detect(variable,"IBD") ~ "IBD",
                             str_detect(variable,"IGT") ~ "IGT",str_detect(variable,"T2D") ~ "T2D")) |>
  mutate(eval = case_when(str_detect(variable,"Accuracy") ~ "Accuracy",str_detect(variable,"F1") ~ "F1",
                          str_detect(variable,"AUROC") ~ "AUROC",str_detect(variable,"AUPRC") ~ "AUPRC")) |>
  mutate(method_eval = paste(method,eval)) |> 
  ggplot(aes(x=value, y=method,fill=method)) +
  geom_rect(aes(ymin=0.5,ymax=2.5,xmin=-Inf,xmax=Inf,fill = "zXGB"),alpha = 0.9,show.legend = FALSE) +
  geom_rect(aes(ymin=4.5,ymax=6.5,xmin=-Inf,xmax=Inf,fill = "zXGB"),alpha = 0.9,show.legend = FALSE) +
  geom_boxplot(outlier.colour="grey75", outlier.size=1,show.legend = FALSE) +
  scale_fill_manual(values=c("#f768a1","#fde0dd","#f16913","#fdd0a2","#3680c0", "#deebf7","#41ae76","#c7e9b4","gray90"))+
  theme_classic()+
  theme(strip.background = element_rect(fill = "ivory"),
        #axis.text.y = element_blank(), 
        #axis.ticks.y=element_blank()
  )+
  scale_y_discrete(limits = rev)+
  labs(y = "")+
  facet_wrap(. ~ eval,nrow = 4) 
dev.off()

