# TE number per peak for each TF

rm(list = ls())
#library(data.table)
library('dplyr')
options(scipen = 999)

files <- list.files('../repeatmasker_track/intsc_rmsk/',pattern = ".bed$")
data <- lapply(files, function(file) { read.table(paste0('../repeatmasker_track/intsc_rmsk/', file), 
                                                stringsAsFactors = FALSE,
                                                sep = '',
                                                header = F,
                                                comment.char = "#")[,4]}) ##colomn4 for repClass, 5 for repFamily


data_cbind <- do.call(cbind, lapply(lapply(data, unlist), `length<-`, max(lengths(data))))

tfname <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))

colnames(data_cbind) <- tfname
#save(data_cbind,file = "te_matrix.Rdata")


#### get TE number
te_num <- vector()
for (i in 1:length(data)) {
  te_num[i] <- data[[i]] %>% length()
}
names(te_num) <- tfname


############=========CTCF RAD21=========######

crfiles <- list.files('../repeatmasker_track/intsc_rmsk/ctcf.rad21/',pattern = ".bed$")
crdata <- lapply(crfiles, function(file) { read.table(paste0('../repeatmasker_track/intsc_rmsk/ctcf.rad21/', file), 
                                                  stringsAsFactors = FALSE,
                                                  sep = '',
                                                  header = F,
                                                  comment.char = "#")[,4]})


crdata_cbind <- do.call(cbind, lapply(lapply(crdata, unlist), `length<-`, max(lengths(crdata))))

crname <- unlist(lapply(X = crfiles, FUN = function(x) {return(strsplit(x, split = "_rmsk")[[1]][1])}))

colnames(crdata_cbind) <- crname


cr_num <- vector()
for (i in 1:length(crdata)) {
  cr_num[i] <- crdata[[i]] %>% length()
}
names(cr_num) <- crname


all.te_num <- c(te_num,cr_num)
##### TF per peak?
#get peak number from hs_getloops
load("../../../hs_getloops/intrachr/1.intrachr_norm/h.all_peak_num.Rdata")


norm_te <- all.te_num/(nm_peaks)
fivenum(norm_te)
# TAFII     GABPA   SMARCA4       YY1     NANOG 
# 0.1261076 0.2318347 0.2998212 0.3433398 0.7384947 



##### distribution plot
library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "none")
)  

p <- ggplot(as.data.frame(norm_te), aes(x = norm_te))
p + geom_histogram(color="darkblue", fill="lightblue",
                   position="identity", alpha=0.5)+
  # geom_vline(aes(xintercept = round(max(norm_te),2)), 
  #            linetype = "dashed", size = 0.6)+
   # geom_text(aes(x = max(norm_te) , y = 15, 
   #               label = tfname[which.max(norm_te)],angle=60
   #               ))+
   # geom_text(aes(x = norm_te[83], y = 15, 
   #              label = tfname[83],angle=60))+
   # geom_text(aes(x = norm_te[108]+5, y = 25, 
   #              label = tfname[108],angle=60))+
   # geom_text(aes(x = norm_te[62], y = 35, 
   #              label = tfname[62],angle=60))+
   labs(title = "TEs Per Peak Distribution",
        x = "TEs Per Peak",
        y = "Count")
ggsave(filename = "hisogram_te.png",path = ".")  

  


#### representive ons, barplot
library(reshape2)
dat <- melt(norm_te)
dat <- cbind(dat,row.names(dat))
colnames(dat) <- c("num","tf")

theme_set(
  theme_minimal() +
    theme(legend.position = "none")
)  


all.name <- c(tfname,crname)
p <- ggplot(data = filter(dat, dat$tf %in% all.name[grep("CTCF", all.name)] |
                            dat$tf %in% all.name[grep("RAD21", all.name)] | 
                            dat$tf %in% all.name[grep("SMC3", all.name)] | 
                            dat$tf %in% all.name[grep("SMC1A", all.name)] | 
                            dat$tf %in% all.name[grep("THAP11", all.name)]),
            aes(x=tf))
p + geom_bar(aes(y = num),stat="identity",fill = "lightblue")+
  geom_text(aes(y = num,label=round(num,3)), vjust=1.3,
            color="black", size=4)+
  coord_flip()+
  labs(x = "TFs",
       y = "Numer of TEs per Peak")
ggsave(filename = "normed_te_number.png",path = ".")  

