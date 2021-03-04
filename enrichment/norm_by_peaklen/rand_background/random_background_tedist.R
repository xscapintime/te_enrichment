# random background TE enrichment
rm(list = ls())
library('dplyr')
require(GenomicRanges)

options(scipen = 999)

files <- list.files('../../repeatmasker_track/intsc_rmsk/rand_backgd/',pattern = ".bed$")
bgname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))


########===== TE class distribution =====########
classdata <- lapply(files, function(file) 
  { read.table(paste0('../../repeatmasker_track/intsc_rmsk/rand_backgd/', file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4]})


### for the ? issue
# for (i in 1:length(classdata)) {
#   classdata[[i]] <- gsub("?","",classdata[[i]],fixed = T)
# }

for (i in 1:length(classdata)) {
  classdata[[i]] <- classdata[[i]][!grepl("?",classdata[[i]],fixed = T)]
}


### all the TE class name
# all.class <- unlist(lapply(classdata,unique)) %>% unique() #if not enough sample, get it from rmsk
load("./cls_order_bycontent.Rdata")


### construct the matrix
cls_mx <- matrix(data=NA, nrow = length(class_order), ncol = length(classdata),
                 byrow = FALSE, dimnames = list(class_order,bgname))

for (cls in class_order) {
  for (j in 1:length(classdata)) {
    
    fq <- as.data.frame(table(classdata[[j]]))
    
    if (fq %>% nrow() < length(class_order)) {
      
      absc <- matrix(nrow = length(setdiff(class_order,fq$Var1)),
                     ncol = 2)
      absc[,1] <-  setdiff(class_order,fq$Var1)
      absc[,2] <- rep(0, length(setdiff(class_order,fq$Var1)))
      colnames(absc) <- c("Var1", "Freq")
      fq <- rbind(fq,absc)
      
    }
    
    cls_mx[cls,j] <- (fq %>% filter(Var1 == cls))[,2]
    
  }
}

### normed TE class
is.matrix(cls_mx)
cls_mx <- apply(cls_mx,2,as.numeric)
row.names(cls_mx) <- class_order

##### normalize the class number with peak length
bg_peaklen <- read.csv("./rand_len_of_peak.txt",header = F,sep = " ")

norm_cls_mx <- t(t(cls_mx)/bg_peaklen[,2])

### mean normed TE class
mean.bg_norm_cls <- rowMeans(norm_cls_mx)
names(mean.bg_norm_cls) <- class_order

save(mean.bg_norm_cls, file = "mean_bg_norm_cls.Rdata")


