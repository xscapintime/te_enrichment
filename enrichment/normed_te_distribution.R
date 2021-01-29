# TE distribution 
rm(list = ls())
library('dplyr')

options(scipen = 999)

#rmsk <- read.table("../repeatmasker_track/rmsk.txt", header = F)

files <- list.files('../repeatmasker_track/intsc_rmsk/',pattern = ".bed$")
tfname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))

crfiles <- list.files('../repeatmasker_track/intsc_rmsk/ctcf.rad21/',pattern = ".bed$")
crname <- unlist(lapply(X = crfiles, FUN = function(x) {return(strsplit(x, split = "_rmsk")[[1]][1])}))

all.name <- c(tfname,crname)

###########======== class dist. ========###########

classdata <- lapply(files, function(file) 
  { read.table(paste0('../repeatmasker_track/intsc_rmsk/', file), 
                                                  stringsAsFactors = FALSE,
                                                  sep = '',
                                                  header = F,
                                                  comment.char = "#")[,4]}) ##column4 for repClass, 5 for repFamily

crclassdata <- lapply(crfiles, function(file) 
{ read.table(paste0('../repeatmasker_track/intsc_rmsk/ctcf.rad21/', file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4]})

all.classdata <- c(classdata,crclassdata)

## for the ? issue
for (i in 1:length(all.classdata)) {
  all.classdata[[i]] <- gsub("?","",all.classdata[[i]],fixed = T)
}

## all the TE class name
all.class <- unlist(lapply(all.classdata,unique)) %>% unique()


## construct matrix
cls_mx <- matrix(data=NA, nrow = length(all.class), ncol = length(all.classdata),
                 byrow = FALSE, dimnames = list(all.class,all.name))

for (cls in all.class) {
  for (j in 1:length(all.classdata)) {
    
    fq <- as.data.frame(table(all.classdata[[j]]))
    
    if (fq %>% nrow() < length(all.class)) {
      absc <- matrix(nrow = length(setdiff(all.class,fq$Var1)),
                     ncol = 2)
      absc[,1] <-  setdiff(all.class,fq$Var1)
      absc[,2] <- rep(0, length(setdiff(all.class,fq$Var1)))
      colnames(absc) <- c("Var1", "Freq")
      fq <- rbind(fq,absc)
    }
    cls_mx[cls,j] <- (fq %>% filter(Var1 == cls))[,2]
  }
}

save(cls_mx, file = "class_matrix.Rdata")

# clsdata_cbind <- do.call(cbind, lapply(lapply(classdata, unlist), `length<-`, max(lengths(classdata))))
# colnames(clsdata_cbind) <- tfname

##### normalize the class number with peak number
load("../../../hs_getloops/intrachr/1.intrachr_norm/h.all_peak_num.Rdata")

is.matrix(cls_mx)
cls_mx <- apply(cls_mx,2,as.numeric)
row.names(cls_mx) <- all.class
norm_cls_mx <- t(t(cls_mx)/nm_peaks)


##### subtract rand background
load("./mean_bg_norm_cls.Rdata")
sbt.norm_cls_mx <- norm_cls_mx - mean.bg_norm_cls[all.class]

save(sbt.norm_cls_mx, file = "all.sbt_norm_cls_mx.Rdata")




