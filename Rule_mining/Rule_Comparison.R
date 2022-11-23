#Rule Comparision

IBD_univ_rule <- fread("univ_rule_py_IBD_supp0.7_uncleaned.csv") 
H_univ_rule <- fread("univ_rule_H_supp0.7_uncleaned.csv") 
IBD_univ_rule$rule_inv <- paste(IBD_univ_rule$V3, IBD_univ_rule$V2)
H_univ_rule$rule_inv <- paste(H_univ_rule$V3, H_univ_rule$V2)
H_IBD_common <- intersect(H_univ_rule$rule,IBD_univ_rule$rule) 
all(intersect(H_univ_rule$rule_inv,IBD_univ_rule$rule) %in% H_IBD_common)
all(H_univ_rule$rule_inv %in% H_univ_rule$rule)
all(IBD_univ_rule$rule_inv %in% IBD_univ_rule$rule)

H_IBD_half_intersect <- H_univ_rule[intersect(which(H_univ_rule$V2 %in% IBD_univ_rule$V2),
                                              which(H_univ_rule$V3 %in% IBD_univ_rule$V3)),]
IBD_H_half_intersect <- IBD_univ_rule[intersect(which(IBD_univ_rule$V2 %in% H_univ_rule$V2),
                                                which(IBD_univ_rule$V3 %in% H_univ_rule$V3)),]


labels_y_H <- levels(as.factor(H_IBD_half_intersect$V3))
breaks_y_H <- seq_along(labels_y_H)

labels_x_H <- levels(as.factor(H_IBD_half_intersect$V2))
breaks_x_H <- seq_along(labels_x_H)

x_H <- as.numeric(factor(H_IBD_half_intersect[H_IBD_half_intersect$rule %in% H_IBD_common,]$V2,
                         levels = c(unique(H_IBD_half_intersect$V2),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V2))))
y_H <- as.numeric(factor(H_IBD_half_intersect[H_IBD_half_intersect$rule %in% H_IBD_common,]$V3,
                         levels = c(unique(H_IBD_half_intersect$V3),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V3))))

p_H_IBD <-  as.data.table(H_IBD_half_intersect) |>
  ggplot(aes(as.integer(factor(V2,levels = c(unique(H_IBD_half_intersect$V2),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V2)))),
             as.integer(factor(V3,levels = c(unique(H_IBD_half_intersect$V3),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V3)))), 
             fill = lift)) +
  geom_tile(aes(fill = lift), color = "black",lwd = 0.05,linetype = 1,show.legend = FALSE) +
  #theme_minimal() +  
  scale_fill_gradient(low = "#d1e5f0", high = "#045a8d") +
  scale_y_continuous(breaks = breaks_y_H, labels = labels_y_H) +
  scale_x_continuous(breaks = breaks_x_H, labels = labels_x_H) +
  theme_bw()+
  xlim(range(as.numeric(factor(c(unique(H_IBD_half_intersect$V2))))))+
  ylim(range(as.numeric(factor(c(unique(H_IBD_half_intersect$V3))))))+
  theme(panel.background = element_rect(fill = 'grey98'),
        panel.grid = element_line(), 
        legend.position="bottom", 
        text = element_text(size = 20),
        axis.text = element_blank(),
        #axis.text.x = element_text(angle = 0, vjust = 0, hjust=0),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour="gray90", size=0.1))+
  labs(title = "", x = "lhs", y = "rhs") + 
  coord_cartesian(expand = FALSE) +
  geom_rect(aes(xmin = c(x_H - offset,rep(0,length(V2)-length(x_H))),
                xmax = c(x_H + offset,rep(0,length(V2)-length(x_H))),
                ymin = c(y_H - offset,rep(0,length(V2)-length(y_H))),
                ymax = c(y_H + offset,rep(0,length(V2)-length(y_H)))),
            fill = "transparent", color = "#d6604d", size = 0.15)

labels_y <- levels(as.factor(IBD_H_half_intersect$V3))
breaks_y <- seq_along(labels_y)

labels_x <- levels(as.factor(IBD_H_half_intersect$V2))
breaks_x <- seq_along(labels_x)

x <- as.numeric(factor(IBD_H_half_intersect[IBD_H_half_intersect$rule %in% H_IBD_common,]$V2,
                       levels = c(unique(H_IBD_half_intersect$V2),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V2))))
y <- as.numeric(factor(IBD_H_half_intersect[IBD_H_half_intersect$rule %in% H_IBD_common,]$V3,
                       levels = c(unique(H_IBD_half_intersect$V3),setdiff(IBD_H_half_intersect$V3,H_IBD_half_intersect$V3))))

offset <- 0.5

p_IBD_H <- as.data.table(IBD_H_half_intersect) |>
  ggplot(aes(as.integer(factor(V2,levels = c(unique(H_IBD_half_intersect$V2),setdiff(V3,H_IBD_half_intersect$V2)))),
             as.integer(factor(V3,levels = c(unique(H_IBD_half_intersect$V3),setdiff(V3,H_IBD_half_intersect$V3)))), fill = lift)) +
  geom_tile(aes(fill = lift), color = "black",lwd = 0.05,linetype = 1,show.legend = FALSE) +
  #theme_minimal() + 
  scale_fill_gradient(low = "#d1e5f0", high = "#045a8d") +
  scale_y_continuous(breaks = breaks_y, labels = labels_y) +
  scale_x_continuous(breaks = breaks_x, labels = labels_x) +
  theme_bw()+
  xlim(range(as.numeric(factor(c(unique(H_IBD_half_intersect$V2))))))+
  ylim(range(as.numeric(factor(c(unique(H_IBD_half_intersect$V3))))))+
  theme(panel.background = element_rect(fill = 'grey98'),
        panel.grid = element_line(), 
        legend.position="bottom", 
        text = element_text(size = 20),
        axis.text = element_blank(),
        #axis.text.x = element_text(angle = 0, vjust = 0, hjust=0),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour="gray90", size=0.1))+
  labs(title = "", x = "lhs", y = "rhs") + 
  coord_cartesian(expand = FALSE) +
  geom_rect(aes(xmin = c(x - offset,rep(0,length(V2)-length(x))),
                xmax = c(x + offset,rep(0,length(V2)-length(x))),
                ymin = c(y - offset,rep(0,length(V2)-length(x))),
                ymax = c(y + offset,rep(0,length(V2)-length(x)))),
            fill = "transparent", color = "#d6604d", size = 0.15)


library(ggpubr)
pdf("Itemset_compare_H_IBD.pdf",width = 15,height = 5)
ggarrange(p_H_IBD, p_IBD_H,
          labels = c("H", "IBD"),
          ncol = 2, nrow = 1,
          font.label=list(color="black",size=20))
dev.off()
