# TE distribution 
rm(list = ls())
library('dplyr')
library(reshape2)

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




##########################################
######### bar chart distribution #########
##########################################
library(ggplot2);
library(RColorBrewer)

theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)
#mycolors <- colorRampPalette(brewer.pal(8,"Set2"))(16)

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
            "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
            "#E34234", "#9A2A2A", "#630330", "#4a3970",
            "#3a58b2", "#607D8B", "#767676", "#36454F")


dat <- melt(sbt.norm_cls_mx) %>% na.omit()
colnames(dat) <- c("te_class","tf","norm")

### fill
p <- ggplot(data = dat %>% filter(norm > 0)
            , aes(x=tf, y=norm, fill=te_class))
p + geom_bar(stat="identity",position="fill")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Class Proportion')+ 
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "class_dist.fill.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "class_dist.fill.pdf", width = 7,height = 12, dpi = 1000)


### stack
q <- ggplot(data = dat
            , aes(x=tf, y=norm, fill=te_class))
q + geom_bar(stat="identity",position="stack")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Class Number')+ 
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "class_dist.stack.png", width = 7,height = 12, dpi = 1000)
ggsave(filename = "class_dist.stack.pdf", width = 7,height = 12, dpi = 1000)


##########################################
############## top and bott ##############
##########################################
top <- read.table("./top.txt")
bott <- read.table("./bott.txt")
tnb <- rbind(top,bott)[,1]

### fill
p2 <- ggplot(data = dat %>% filter(norm > 0 & (tf %in% tnb) == T)
         , aes(x=tf, y=norm, fill=te_class))
p2 + geom_bar(stat="identity",position="fill")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Class Proportion')+ 
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "tnb-class_dist.fill.png", width = 7,height = 9, dpi = 1000)
ggsave(filename = "tnb-class_dist.fill.pdf", width = 7,height = 9, dpi = 1000)


### stack
q2 <- ggplot(data = dat %>% filter(tf %in% tnb == T)
            , aes(x=tf, y=norm, fill=te_class))
q2 + geom_bar(stat="identity",position="stack")+
  labs(title = "TE Class Distribution",
       x = NULL,y = 'Normed TE Class Number')+ 
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))
ggsave(filename = "tnb-class_dist.stack.png", width = 7,height = 9, dpi = 1000)
ggsave(filename = "tnb-class_dist.stack.pdf", width = 7,height = 9, dpi = 1000)



