### stack bar plot for looping zscore top22(+3 RAD21 +4 CTCF) and bott29
rm(list = ls())
library(dplyr)
library(reshape2)

load("te_matrix.Rdata")
bott <- read.table("bott29.txt")[,1]
top <- read.table("top22.txt")[,1]


top_list <- list()
for (i in 1:length(top)) {
  top_list[[i]] <- colnames(data_cbind)[grep(paste0(top[i],"_"), 
                                             colnames(data_cbind))]
}
top_list <- unlist(top_list)


bott_list <- list()
for (i in 1:length(bott)) {
  bott_list[[i]] <- colnames(data_cbind)[grep(paste0(bott[i],"_"), 
                                             colnames(data_cbind))]
}
bott_list <- unlist(bott_list)



#### stack bar plot
library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "none")
)

top_dat <- na.omit(melt(data_cbind[,top_list]))
colnames(top_dat)[2:3] <- c("sample","te")
save(top_dat,file = "te_enrich_topzscore-tf.Rdata")
 
topp <- ggplot(data = top_dat, aes(x=sample, y=te,fill=te))
topp + geom_bar(stat="identity",position="stack") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


