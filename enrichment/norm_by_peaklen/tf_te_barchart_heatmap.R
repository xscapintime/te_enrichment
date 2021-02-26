#### peak number normed matrix vis. ####
rm(list = ls())
library(dplyr)
library(reshape2)

options(scipen = 999)

### load the matrix
load("./nromedpkn_class_matrix_sorted.Rdata")

### transform to long data format
dat <- melt(normpkn_cls_mx)
colnames(dat) <- c("class", "tf", "norm")

## looping zscore
lp_zsc <- read.csv("../.././../../hs_getloops/intrachr/3.plot/exless_h.all.tf_des-meanz.v3.csv", header = T)
# group
lp_zsc$group <- ifelse(lp_zsc$allmean_z >= 0.8, "high",
                       ifelse(lp_zsc$allmean_z <= -0.8, "low",
                              "mid"))
# merge the group
dat <- merge(dat, lp_zsc$group)

## label the tf with colors
dat$tf <- factor(dat$tf, levels = rev(colnames(normpkn_cls_mx))) #reverse the tf order

labelcol <-ifelse(lp_zsc$group == "high","#CD5C5C",   #'high' = 'IndianRed'
                  ifelse(lp_zsc$group == "mid", "#DAA520", #'mid' = 'Goldenrod'
                         "#4682B4"))  #'low' = 'SteelBlue'

## barchart to show distribution
library(ggplot2)

theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)

library(ggsci)
simpsons <- pal_simpsons("springfield")(16) 

# fill
p <- ggplot(data = dat,
            aes(x=tf, y=norm, fill=class))
p + geom_bar(stat="identity",position="fill")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Proportion')+ 
  coord_flip()+
  scale_fill_manual(values = simpsons)+
  guides(fill = guide_legend(title="TE Class", reverse = F))+
  theme(axis.text.x = element_text(size = 8,),
        axis.text.y = element_text(size = 7, colour = rev(labelcol), face = "bold"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "../../figure/norm_resize/tf_class_bar_fill.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "../../figure/norm_resize/tf_class_bar_fill.pdf", width = 7,height = 12, dpi = 1000)

# stack
p <- ggplot(data = dat,
            aes(x=tf, y=norm, fill=class))
p + geom_bar(stat="identity",position="stack")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Number')+ 
  coord_flip()+
  scale_fill_manual(values = simpsons)+
  guides(fill = guide_legend(title="TE Class", reverse = F))+
  theme(axis.text.x = element_text(size = 8,),
        axis.text.y = element_text(size = 7, colour = rev(labelcol), face = "bold"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "../../figure/norm_resize/tf_class_bar_stack.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "../../figure/norm_resize/tf_class_bar_stack.pdf", width = 7,height = 12, dpi = 1000)


## histogram

theme_set(
  theme_classic() +
    theme(legend.position = "none")
)  

p <- ggplot(dat, aes(x = norm))
p + geom_histogram(aes(y = stat(density),colour = "black"), 
                   fill="white",binwidth = 0.003) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  scale_color_manual(values = c("#868686FF")) +
  labs(title = "Normalized TE Num Histogram",
       subtitle = "Human",
       caption = "binwidth = 0.003") +
  ylab('Dens.') +
  xlab('Normalized TE number')
  
ggsave(filename = "Normalized_te_number_hisogram.png",path = "../../figure/norm_resize")  




########### next time start form here #############
## set the break, then heatmap

range(norm_cls_mx*200) #resized to 200 bp, back to the peak number
fivenum(norm_cls_mx*200)
bk <- seq(min(norm_cls_mx*200), max(norm_cls_mx*200), by = 0.003)

#library("RColorBrewer")
library(viridis)
#col <- colorRampPalette(brewer.pal(11, "Spectral"))(length(bk)) %>% rev()

# default heatmap
pheatmap::pheatmap(t(norm_cls_mx*200), scale = "none", 
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = bk,
                   color = viridis(length(bk)),
                   #border_color = "white",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap default",
                   cellwidth = 12,
                   cellheight = 12,
                   filename = "../../figure/norm_resize/default_te_heatmap.pdf")

pheatmap::pheatmap(t(norm_cls_mx*200), scale = "row", 
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = bk,
                   color = viridis(length(bk)),
                   #border_color = "white",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap scaled by row",
                   cellwidth = 12,
                   cellheight = 12,
                   filename = "../../figure/norm_resize/rowscale_te_heatmap.pdf")

pheatmap::pheatmap(t(norm_cls_mx*200), scale = "column", 
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = bk,
                   color = viridis(length(bk)),
                   #border_color = "white",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap scaled by column",
                   cellwidth = 12,
                   cellheight = 12,
                   filename = "../../figure/norm_resize/columnscale_te_heatmap.pdf")



