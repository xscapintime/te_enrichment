rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

## load peak files
peakfiles <- list.files('../../../hs_getloops/data/peak/',pattern = ".bed$")
peak <- lapply(peakfiles, function(file) { import.bed(paste0('../../../hs_getloops/data/peak/', file))})
peakname <- unlist(lapply(X = peakfiles, FUN = function(x) {return(strsplit(x, split = "_peaks")[[1]][1])}))
names(peak) <- peakname

peaks_resized = lapply(peak, function(p) {resize(p,width = 200, fix = 'center')})

## export resized peak files
for (i in 1:length(peaks_resized)) {
  export(peaks_resized[[i]], paste0("resized_peaks/",peakname[i],"_resized.bed"))
}



### for CTCF and RAD21
crpeakfiles <- list.files('../../../hs_getloops/data/ctcf.rad21/peak/',pattern = "bed$")
crpeak <- lapply(crpeakfiles, function(file) { import.bed(paste0('../../../hs_getloops/data/ctcf.rad21/peak/', file))})
crpeakname <- unlist(lapply(X = crpeakfiles, FUN = function(x) {return(strsplit(x, split = "_peaks")[[1]][1])}))
names(crpeak) <- crpeakname

crpeaks_resized = lapply(crpeak, function(p) {resize(p,width = 200, fix = 'center')})


### export resized peak files
for (j in 1:length(crpeaks_resized)) {
  export(crpeaks_resized[[j]], paste0("resized_peaks/ctcf.rad21/",crpeakname[j],"_resized.bed"))
}
