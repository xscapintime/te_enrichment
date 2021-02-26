#### try quantile normalization ####
# try quantile normalizaiton
#BiocManager::install("preprocessCore")
library(preprocessCore)

qn_mx <- normalize.quantiles(t(norm_cls_mx*200))
dimnames(qn_mx) <- list(colnames(norm_cls_mx),row.names(norm_cls_mx))



## quntail breaks?
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(norm_cls_mx*200, n =  10)



## first barchat to see the difference

## 
pheatmap::pheatmap(qn_mx, scale = "none", 
                   cluster_cols = T,
                   cluster_rows = F,
                   #breaks = bk,
                   color = viridis(length(bk)),
                   #border_color = "white",
                   clustering_distance_rows = "correlation",
                   main = "pheatmap quantile normalized",
                   cellwidth = 12,
                   cellheight = 12,
                   filename = "../../figure/norm_resize/qn_te_heatmap.pdf")



