###### te number normed by total peak length(unified) ######
rm(list = ls())
library(dplyr)

options(scipen = 999)

files <- list.files('../../repeatmasker_track/resized_intsc/',pattern = ".bed$")
tfname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))

crfiles <- list.files('../../repeatmasker_track/resized_intsc/ctcf.rad21/',pattern = ".bed$")
crname <- unlist(lapply(X = crfiles, FUN = function(x) {return(strsplit(x, split = "_resized")[[1]][1])}))

all.name <- c(tfname,crname)

###########======== load intersected data ========###########

classdata <- lapply(files, function(file) 
{ read.table(paste0('../../repeatmasker_track/resized_intsc/', file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4]}) ##column4 for repClass, 5 for repFamily

crclassdata <- lapply(crfiles, function(file) 
{ read.table(paste0('../../repeatmasker_track/resized_intsc/ctcf.rad21/', file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4]})

all.classdata <- c(classdata,crclassdata)

## all the TE class name (decreasing order by content)
load("te_cls_count.Rdata")
class_order <- te_cls_count$V4

## construct matrix
cls_mx <- matrix(data=NA, nrow = length(class_order), ncol = length(all.classdata),
                 byrow = FALSE, dimnames = list(class_order,all.name))

for (cls in class_order) {
  for (j in 1:length(all.classdata)) {
    
    fq <- as.data.frame(table(all.classdata[[j]]))
    
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

## tf order(looping zscore)
lp_zsc <- read.csv("../.././../../hs_getloops/intrachr/3.plot/exless_h.all.tf_des-meanz.v3.csv", header = T)

cls_mx <- cls_mx[,lp_zsc[,1]]
save(cls_mx, file = "class_matrix_sorted.Rdata")


### normed by unified peak length 
#since the peak length is unified(all 200 bp), maybe normed by peak number?
unied_peaklen <- read.csv("./len_of_resizedpeak.txt", sep = " ", header = F, row.names = 1)
peak_len <- unied_peaklen[lp_zsc[,1],][,1]
names(peak_len) <- lp_zsc[,1]

is.matrix(cls_mx)
cls_mx <- apply(cls_mx,2,as.numeric)
row.names(cls_mx) <- class_order

norm_cls_mx <- t(t(cls_mx)/peak_len)

### normed by peak number
normpkn_cls_mx <- norm_cls_mx*200
save(normpkn_cls_mx, file = "nromedpkn_class_matrix_sorted.Rdata") 

