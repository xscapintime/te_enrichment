# TE(tf)/TF(g)/peaks(tf) - TE(rdg)/TF(g)/peaks(rdg) = [TE(tf)/p(tf) - TE(rdg)/p(rdg)]TE(genome)
rm(list = ls())
library(dplyr)

## class
### TE in the genome
gn.te_cls_mx <- readRDS("../teclass_mx_bychr.rds")

### normed TE (divide peaks and subtract rand bg)
load("../all.sbt_norm_cls_mx.Rdata")

### [TE(tf)/p(tf) - TE(rdg)/p(rdg)]TE(genome)
all.class <- rownames(gn.te_cls_mx)
gn.te_cls_mx <- apply(gn.te_cls_mx, 2, as.numeric)
row.names(gn.te_cls_mx) <- all.class

gn.te_cls_num <-  apply(gn.te_cls_mx, 1, sum)

cls_normed_by_gtenum <- sbt.norm_cls_mx[all.class,]/gn.te_cls_num
#sbt.norm_cls_mx[all.class,][,2]/gn.te_cls_num

############################################
####### TE class distribution for TF #######
############################################
library(ggplot2)
library(reshape2)
library(RColorBrewer)

dat <- melt(cls_normed_by_gtenum)
colnames(dat) <- c("class", "tf", "normed")

theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
) 

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
            "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
            "#E34234", "#9A2A2A", "#630330", "#342D7E",
            "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

p <- ggplot(dat %>% filter(class %in% c("Unknown","tRNA") == F &
                             normed > 0), 
            aes(x=tf, y=normed, fill=class))
p + geom_col(position="fill")+
  labs(title = "TE Class Distribution in TF Binding Sites",
       x = NULL,y = 'log10 TE Class Proportion')+
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 8),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()
