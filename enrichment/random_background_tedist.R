# random background TE enrichment
rm(list = ls())
library('dplyr')
require(GenomicRanges)

options(scipen = 999)

files <- list.files('../repeatmasker_track/intsc_rmsk/rand_backgd/',pattern = ".bed$")
bgname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))

########===== TE number per peak =====########
classdata <- lapply(files, function(file) 
  { read.table(paste0('../repeatmasker_track/intsc_rmsk/rand_backgd/', file), 
                                                       stringsAsFactors = FALSE,
                                                       sep = '',
                                                       header = F,
                                                       comment.char = "#")[,4]})

te_num <- vector()
for (i in 1:length(classdata)) {
  te_num[i] <- classdata[[i]] %>% length()
}
names(te_num) <- bgname

### get peak number from hs_getloops/
load("../../../hs_getloops/intrachr/2.intrachr_sbtr/h.presumed_peaks.Rdata")
peak_num <- lengths(peak)

### normalized TE number
norm_te <- te_num/peak_num
### mean normed background TE number 
m.norm_te <- mean(norm_te)
print(m.norm_te)


########===== TE class distribution =====########
classdata <- lapply(files, function(file) 
  { read.table(paste0('../repeatmasker_track/intsc_rmsk/rand_backgd/', file), 
             stringsAsFactors = FALSE,
             sep = '',
             header = F,
             comment.char = "#")[,4]})


### for the ? issue
for (i in 1:length(classdata)) {
  classdata[[i]] <- gsub("?","",classdata[[i]],fixed = T)
}

### all the TE class name
all.class <- unlist(lapply(classdata,unique)) %>% unique() #if not enough sample, get it from rmsk

### construct the matrix
cls_mx <- matrix(data=NA, nrow = length(all.class), ncol = length(classdata),
                 byrow = FALSE, dimnames = list(all.class,bgname))

for (cls in all.class) {
  for (j in 1:length(classdata)) {
    
    fq <- as.data.frame(table(classdata[[j]]))
    
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

### normed TE class
is.matrix(cls_mx)
cls_mx <- apply(cls_mx,2,as.numeric)
row.names(cls_mx) <- all.class
norm_cls_mx <- cls_mx/peak_num

### mean normed TE class
mean.bg_norm_cls <- rowMeans(norm_cls_mx)
names(mean.bg_norm_cls) <- all.class

save(mean.bg_norm_cls, file = "mean_bg_norm_cls.Rdata")


