###########======== family dist. ========###########
rm(list = ls())

files <- list.files('../repeatmasker_track/intsc_rmsk/',pattern = ".bed$")
tfname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))


fmldata <- lapply(files, function(file) { read.table(paste0('../repeatmasker_track/intsc_rmsk/', file), 
                                                     stringsAsFactors = FALSE,
                                                     sep = '',
                                                     header = F,
                                                     comment.char = "#")[,5]}) ##column4 for repClass, 5 for repFamily

## for the ? issue
for (i in 1:length(fmldata)) {
  fmldata[[i]] <- gsub("?","",fmldata[[i]],fixed = T)
}

## all the TE family name
all.fml <- unlist(lapply(fmldata,unique)) %>% unique()

## construct matrix
fml_mx <- matrix(data=NA, nrow = length(all.fml), ncol = length(fmldata),
                 byrow = FALSE, dimnames = list(all.fml,tfname))

for (cls in all.fml) {
  for (j in 1:length(fmldata)) {
    
    fq <- as.data.frame(table(fmldata[[j]]))
    
    if (fq %>% nrow() < length(all.fml)) {
      absc <- matrix(nrow = length(setdiff(all.fml,fq$Var1)),
                     ncol = 2)
      absc[,1] <-  setdiff(all.fml,fq$Var1)
      absc[,2] <- rep(0, length(setdiff(all.fml,fq$Var1)))
      colnames(absc) <- c("Var1", "Freq")
      fq <- rbind(fq,absc)
    }
    fml_mx[cls,j] <- (fq %>% filter(Var1 == cls))[,2]
  }
}

save(fml_mx, file = "family_matrix.Rdata")


#fmldata_cbind <- do.call(cbind, lapply(lapply(fmldata, unlist), `length<-`, max(lengths(fmldata))))
#colnames(fmldata_cbind) <- tfname
