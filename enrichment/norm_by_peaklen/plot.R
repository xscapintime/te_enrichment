# plot

rm(list = ls())
library('dplyr')

options(scipen = 999)

# tf looping order
tf_order <- read.csv("../../../../hs_getloops/intrachr/3.plot/exless_h.all.tf_des-meanz.v2.csv")

# barplot
## norm by peak length
load("./cls_order_bycontent.Rdata")
load("./norm_cls_mx.Rdata")

library(reshape2)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)

# mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
#             "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
#             "#E34234", "#9A2A2A", "#630330", "#4a3970",
#             "#3a58b2", "#607D8B", "#767676", "#36454F") %>% rev()
# 
# mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
#             "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
#             "#E34234", "#9A2A2A", "#630330", "#342D7E",
#             "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

library(ggsci)
simpsons <- pal_simpsons("springfield")(16) %>% rev()


dat <- melt(norm_cls_mx[rev(class_order),]) %>% na.omit()
colnames(dat) <- c("te_class","tf","norm")
p <- ggplot(data = dat,
            aes(x=tf, y=norm, fill=te_class))
p + geom_bar(stat="identity",position="stack")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Class Proportion')+ 
  coord_flip()+
  scale_fill_manual(values = simpsons)+
  guides(fill = guide_legend(title="TE Class", reverse = T))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))


