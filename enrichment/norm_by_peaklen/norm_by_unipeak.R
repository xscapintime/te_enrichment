###### te number normed by total peak length(unified) ######
rm(list = ls())
library(dplyr)

options(scipen = 999)

getwd()
setwd("./human/enrichment/norm_by_peaklen")

files <- list.files("../../repeatmasker_track/resized_intsc/", pattern = ".bed$")
tfname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))

crfiles <- list.files("../../repeatmasker_track/resized_intsc/ctcf.rad21/", pattern = ".bed$")
crname <- unlist(lapply(X = crfiles, FUN = function(x) {return(strsplit(x, split = "_resized")[[1]][1])}))

all.name <- c(tfname,crname)

## tf order(looping zscore)
lp_zsc <- read.csv("../.././../../hs_getloops/intrachr/3.plot/exless_h.all.tf_des-meanz.v3.csv", header = T)

#### TE class
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
  for (j in seq_along(all.classdata)) {
    
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

cls_mx <- cls_mx[,lp_zsc[,1]]
save(cls_mx, file = "class_matrix_sorted.Rdata")

### normed by unified peak length
#since the peak length is unified(all 200 bp), maybe normed by peak number?
unied_peaklen <- read.csv("./len_of_resizedpeak.txt",
                          sep = " ", header = F, row.names = 1)
peak_len <- unied_peaklen[lp_zsc[,1],][,1]
names(peak_len) <- lp_zsc[,1]

is.matrix(cls_mx)
cls_mx <- apply(cls_mx, 2, as.numeric)
row.names(cls_mx) <- class_order

norm_cls_mx <- t(t(cls_mx) / peak_len)

### normed by peak number
normpkn_cls_mx <- norm_cls_mx*200
save(normpkn_cls_mx, file = "nromedpkn_class_matrix_sorted.Rdata") 

#### TE family
###########======== load intersected data ========###########

fmldata <- lapply(files, function(file)
{ read.table(paste0("../../repeatmasker_track/resized_intsc/", file),
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4:5]}) ##column4 for repClass, 5 for repFamily

crfmldata <- lapply(crfiles, function(file)
{ read.table(paste0("../../repeatmasker_track/resized_intsc/ctcf.rad21/", file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4:5]})

all.fmldata <- c(fmldata, crfmldata)

## merge the class and family name together
for (i in seq_along(all.fmldata)) {
  all.fmldata[[i]] <- all.fmldata[[i]] %>% mutate(fml = paste(V4, V5))
  all.fmldata[[i]] <- all.fmldata[[i]]$fml
}

## all the TE family name (decreasing order by content)
load("te_fml_count.Rdata")
fml_order <- te_fml_count$V5

## construct matrix
fml_mx <- matrix(data = NA, nrow = length(fml_order), ncol = length(all.fmldata),
                 byrow = FALSE, dimnames = list(fml_order, all.name))

for (fml in fml_order) {
  for (j in seq_along(all.fmldata)) {
    
    fq <- as.data.frame(table(all.fmldata[[j]]))
    
    if (fq %>% nrow() < length(fml_order)) {
      absc <- matrix(nrow = length(setdiff(fml_order, fq$Var1)),
                     ncol = 2)
      absc[,1] <-  setdiff(fml_order, fq$Var1)
      absc[,2] <- rep(0, length(setdiff(fml_order, fq$Var1)))
      colnames(absc) <- c("Var1", "Freq")
      fq <- rbind(fq, absc)
    }
    fml_mx[fml,j] <- (fq %>% filter(Var1 == fml))[,2]
  }
}

fml_mx <- fml_mx[,lp_zsc[,1]]
save(fml_mx, file = "family_matrix_sorted.Rdata")

### normed by unified peak length
#since the peak length is unified(all 200 bp), maybe normed by peak number?
unied_peaklen <- read.csv("./len_of_resizedpeak.txt",
                          sep = " ", header = F, row.names = 1)
peak_len <- unied_peaklen[lp_zsc[,1],][,1]
names(peak_len) <- lp_zsc[,1]

is.matrix(fml_mx)
fml_mx <- apply(fml_mx, 2, as.numeric)
row.names(fml_mx) <- fml_order

norm_fml_mx <- t(t(fml_mx) / peak_len)

### normed by peak number
normpkn_fml_mx <- norm_fml_mx * 200
save(normpkn_fml_mx, file = "nromedpkn_family_matrix_sorted.Rdata")
