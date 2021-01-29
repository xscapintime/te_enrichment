# TE distribution in TF binding sites by class/family & plots

rm(list = ls())
library(dplyr)
library(reshape2)
load("all.sbt_norm_cls_mx.Rdata")

##########################################
######### bar chart distribution #########
##########################################
library(ggplot2);
library(RColorBrewer)

theme_set(
	    theme_minimal() +
		        theme(legend.position = "right",
			                legend.key.size = unit(12, "pt"))
	    )
#mycolors <- colorRampPalette(brewer.pal(8,"Set2"))(16)

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
	                "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
			            "#E34234", "#9A2A2A", "#630330", "#4a3970",
			            "#3a58b2", "#607D8B", "#767676", "#36454F")


dat <- melt(sbt.norm_cls_mx) %>% na.omit()
colnames(dat) <- c("te_class","tf","norm")

### fill
p <- ggplot(data = dat %>% filter(norm > 0)
	                , aes(x=tf, y=norm, fill=te_class))
p + geom_bar(stat="identity",position="fill")+
	  labs(title = "TE Class Distribution",
	              x = NULL,y = 'Normed TE Class Proportion')+ 
  coord_flip()+
    scale_fill_manual(values = mycols)+
      guides(fill = guide_legend(title="TE Class"))+
        theme(axis.text.x = element_text(size = 8),
	              axis.text.y = element_text(size = 7),
		              legend.title=element_text(size=7),
		              legend.text=element_text(size=6))
ggsave(filename = "class_dist.fill.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "class_dist.fill.pdf", width = 7,height = 12, dpi = 1000)


### stack
q <- ggplot(data = dat
                 , aes(x=tf, y=norm, fill=te_class))
q + geom_bar(stat="identity",position="stack")+
	    labs(title = "TE Class Distribution",
		        x = NULL,y = 'Normed TE Class Number')+ 
  coord_flip()+
    scale_fill_manual(values = mycols)+
      guides(fill = guide_legend(title="TE Class"))+
        theme(axis.text.x = element_text(size = 8),
	              axis.text.y = element_text(size = 7),
		              legend.title=element_text(size=7),
		              legend.text=element_text(size=6))
ggsave(filename = "class_dist.stack.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "class_dist.stack.pdf", width = 7,height = 12, dpi = 1000)


##########################################
############## top and bott ##############
##########################################
top <- read.table("./top.txt")
bott <- read.table("./bott.txt")
tnb <- rbind(top,bott)[,1]

### fill
p2 <- ggplot(data = dat %>% filter(norm > 0 & (tf %in% tnb) == T)
               , aes(x=tf, y=norm, fill=te_class))
p2 + geom_bar(stat="identity",position="fill")+
    labs(title = "TE Class Distribution",
	        x = NULL,y = 'Normed TE Class Proportion')+ 
  coord_flip()+
    scale_fill_manual(values = mycols)+
      guides(fill = guide_legend(title="TE Class"))+
        theme(axis.text.x = element_text(size = 8),
	              axis.text.y = element_text(size = 7.5),
		              legend.title=element_text(size=7),
		              legend.text=element_text(size=6))
ggsave(filename = "tnb-class_dist.fill.png", width = 7,height = 9, dpi = 1000)
ggsave(filename = "tnb-class_dist.fill.pdf", width = 7,height = 9, dpi = 1000)


### stack
q2 <- ggplot(data = dat %>% filter(tf %in% tnb == T)
                   , aes(x=tf, y=norm, fill=te_class))
q2 + geom_bar(stat="identity",position="stack")+
	    labs(title = "TE Class Distribution",
		        x = NULL,y = 'Normed TE Class Number')+ 
  coord_flip()+
    scale_fill_manual(values = mycols)+
      guides(fill = guide_legend(title="TE Class"))+
        theme(axis.text.x = element_text(size = 8),
	              axis.text.y = element_text(size = 7.5),
		              legend.title=element_text(size=7),
		              legend.text=element_text(size=6))
ggsave(filename = "tnb-class_dist.stack.png", width = 7,height = 9, dpi = 1000)
ggsave(filename = "tnb-class_dist.stack.pdf", width = 7,height = 9, dpi = 1000)




