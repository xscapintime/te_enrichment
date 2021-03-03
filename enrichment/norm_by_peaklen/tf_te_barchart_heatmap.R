#### peak number normed matrix vis. ####
rm(list = ls())
getwd()
setwd("./human/enrichment/norm_by_peaklen")
getwd()

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
# merge the group (no need to)
#dat <- merge(dat, lp_zsc$group)

## label the tf with colors
dat$tf <- factor(dat$tf,
                 levels = rev(colnames(normpkn_cls_mx)))  #reverse the tf order

labelcol <- ifelse(lp_zsc$group == "high", "#CD5C5C",   #'high' = 'IndianRed'
                  ifelse(lp_zsc$group == "mid", "#DAA520", #'mid' = 'Goldenrod'
                         "#4682B4"))  #'low' = 'SteelBlue'

### default matrix
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
            aes(x = tf, y = norm, fill = class))
p + geom_bar(stat = "identity", position = "fill") +
  labs(title = "TE Class Distribution",
       x = NULL,y = "Normed TE Proportion") +
  coord_flip() +
  scale_fill_manual(values = simpsons) +
  guides(fill = guide_legend(title="TE Class", reverse = F)) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7, colour = rev(labelcol), face = "bold"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))
ggsave(filename = "../../figure/norm_resize/tf_class_bar_fill.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "../../figure/norm_resize/tf_class_bar_fill.pdf", width = 7,height = 12, dpi = 1000)

# stack
p <- ggplot(data = dat,
            aes(x = tf, y = norm, fill = class))
p + geom_bar(stat = "identity", position = "stack") +
  labs(title = "TE Class Distribution",
       x = NULL, y = "Normed TE Number") +
  coord_flip() +
  scale_fill_manual(values = simpsons) +
  guides(fill = guide_legend(title = "TE Class", reverse = F)) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7, colour = rev(labelcol), face = "bold"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))
ggsave(filename = "../../figure/norm_resize/tf_class_bar_stack.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "../../figure/norm_resize/tf_class_bar_stack.pdf", width = 7,height = 12, dpi = 1000)


## histogram
theme_set(
  theme_classic() +
    theme(legend.position = "none")
)

p <- ggplot(dat, aes(x = norm))
p + geom_histogram(aes(y = stat(density), colour = "black"),
                   fill = "white", bins = 100) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  scale_color_manual(values = c("#868686FF")) +
  labs(title = "Normalized TE Num Histogram",
       subtitle = "Human",
       caption = "bins = 100") +
  ylab("Dens.") +
  xlab("Normalized TE number")
ggsave(filename = "Normalized_te_number_hisogram.png", path = "../../figure/norm_resize")  


## set the break, then heatmap
range(normpkn_cls_mx) #resized to 200 bp, back to the peak number
quantile(normpkn_cls_mx)
bk <- seq(min(normpkn_cls_mx), max(normpkn_cls_mx), length.out = 100) #uniform breaks

# Heatmap!
#library("RColorBrewer")
library(viridis)
#col <- colorRampPalette(brewer.pal(11, "Spectral"))(length(bk)) %>% rev()

## group of the tf
anno_tf <- data.frame(lp_zsc$group, row.names = lp_zsc$X)
colnames(anno_tf) <- "looping score"

## anno color
looping <- list(`looping score` = c(high = "#CD5C5C",
                                    mid = "#DAA520",
                                    low = "#4682B4"))

# default heatmap
pheatmap::pheatmap(t(normpkn_cls_mx), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = bk,
                   color = viridis(length(bk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap default",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/default_te_heatmap_unibk100.png")

# tf cluster
pheatmap::pheatmap(t(normpkn_cls_mx), scale = "none",
                   cluster_cols = T,
                   cluster_rows = T,
                   breaks = bk,
                   color = viridis(length(bk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap default",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/default_te_heatmap_cluster_unibk100.png")

## try quantile breaks
#https://slowkow.com/notes/pheatmap-tutorial/#quantile-breaks
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
qbk <- quantile_breaks(normpkn_cls_mx, n = 100)

# quantile breaks heatmap
pheatmap::pheatmap(t(normpkn_cls_mx), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = qbk,
                   color = viridis(length(qbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap default quantile breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/default_te_heatmap_qbk100.png")

# tf cluster
pheatmap::pheatmap(t(normpkn_cls_mx), scale = "none",
                   cluster_cols = T,
                   cluster_rows = T,
                   breaks = qbk,
                   color = viridis(length(qbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap default quantile breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/default_te_heatmap_cluster_qbk100.png")


### z-score normalization
zn_tecls <- t(scale(t(normpkn_cls_mx))) #(n-row.mean)/row.sd

## long data format for ggplot2 histogram
zdat <- melt(zn_tecls)
colnames(zdat) <- c("class", "tf", "norm")

## histogram
theme_set(
  theme_classic() +
    theme(legend.position = "none")
)

p <- ggplot(zdat, aes(x = norm))
p + geom_histogram(aes(y = stat(density), colour = "black"),
                   fill = "white", bins = 50) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  scale_color_manual(values = c("#868686FF")) +
  labs(title = "Z-Normalized TE Num Histogram",
       subtitle = "Human",
       caption = "bins = 50") +
  ylab("Dens.") +
  xlab("Normalized TE number")
ggsave(filename = "zs-Normalized_te_number_hisogram.png", path = "../../figure/norm_resize")  

## breaks
range(zn_tecls)
quantile(zn_tecls)
zbk <- seq(min(zn_tecls), max(zn_tecls), length.out = 100) #uniform breaks
qzbk <- quantile_breaks(zn_tecls, n = 100) #quantile breaks

# uniform breaks
pheatmap::pheatmap(t(zn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = zbk,
                   color = viridis(length(zbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap z-score normalization uniform breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/zn_te_heatmap_unibk.png")

# uni breaks TFs clustered
pheatmap::pheatmap(t(zn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = T,
                   breaks = zbk,
                   color = viridis(length(zbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap z-score normalization uniform breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/zn_te_heatmap_cluster_unibk.png")

# quantile breaks
pheatmap::pheatmap(t(zn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = qzbk,
                   color = viridis(length(qzbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap z-score normalization quantile breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/zn_te_heatmap_qbk.png")

# qb tf cluster
pheatmap::pheatmap(t(zn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = T,
                   breaks = qzbk,
                   color = viridis(length(qzbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap z-score normalization quantile breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/zn_te_heatmap_cluster_qbk.png")


### quantile normalization
library(preprocessCore)
qn_tecls <- normalize.quantiles(normpkn_cls_mx)
dimnames(qn_tecls) <- list(row.names(zn_tecls), colnames(zn_tecls))

## long data format for ggplot2 histogram
qdat <- melt(qn_tecls)
colnames(qdat) <- c("class", "tf", "norm")

## histogram
theme_set(
  theme_classic() +
    theme(legend.position = "none")
)

p <- ggplot(qdat, aes(x = norm))
p + geom_histogram(aes(y = stat(density), colour = "black"),
                   fill = "white", bins = 18) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  scale_color_manual(values = c("#868686FF")) +
  labs(title = "Quantile Normalized TE Num Histogram",
       subtitle = "Human",
       caption = "bins = 18") +
  ylab("Dens.") +
  xlab("Quantile Normalized TE number")
ggsave(filename = "quantile-Normalized_te_number_hisogram.png", path = "../../figure/norm_resize")

## breaks
range(qn_tecls)
quantile(qn_tecls)
uni_qnbk <- seq(min(qn_tecls), max(qn_tecls), length.out = 18) #uniform breaks
q_qnbk <- quantile_breaks(qn_tecls, n = 18) #quantile breaks

# uniform breaks
pheatmap::pheatmap(t(qn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = uni_qnbk,
                   color = viridis(length(uni_qnbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap quantile normalization uniform breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/qn_te_heatmap_unibk.png")

# quantile breaks
pheatmap::pheatmap(t(qn_tecls), scale = "none",
                   cluster_cols = T,
                   cluster_rows = F,
                   breaks = q_qnbk,
                   color = viridis(length(q_qnbk)),
                   border_color = "grey20",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap quantile normalization quantile breaks",
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = anno_tf,
                   annotation_colors = looping,
                   filename = "../../figure/norm_resize/qn_te_heatmap_qbk.png")

#---
#Seems quantile normalizaiton doesn't work so well
