# TE distribution in human genome
rm(list = ls())
library(dplyr)
library(reshape2)

options(scipen = 999)

## load rmsk table
rmsk <- read.table("../repeatmasker_track/del.qmark_rmsk.txt", header = F)
genome_te <- cbind(rmsk$V6, rmsk$V12, paste(rmsk$V12,rmsk$V13)) %>% as.data.frame()


## include the alt contig 
### way too slow
# for (i in grep("_",genome_te$chr)) {
#   
#     genome_te$chr[i] <- strsplit(genome_te$chr[i], "_")[[1]][1]
#   
# }

chr_name <- read.table("../repeatmasker_track/chr_name.txt")
genome_te <- cbind(genome_te, chr_name)
colnames(genome_te) <- c("chr", "class", "family","chr_inc")

table(genome_te$class) #no ? occure
table(genome_te$family) #no ? occure

all.class <- unique(genome_te$class) %>% sort()
all.fml <- unique(genome_te$family) %>% sort()

save(genome_te,file = "genome_te.Rdata")

##########################################
######## construct matrix | class ########
##########################################
cls_mx <- matrix(data=NA, nrow = length(all.class), ncol = length(genome_te$chr_inc %>% unique()),
                 byrow = FALSE, dimnames = list(all.class,genome_te$chr_inc %>% unique()))

for (cls in all.class) {
  for (j in 1:length(genome_te$chr_inc %>% unique())) {
    
    fq <- as.data.frame(table(filter(genome_te,
                               chr_inc==unique(genome_te$chr_inc)[j])[,"class"]))
    
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


#################=============== plot ================#########################



##########################################
####### class pie chart for genome #######
##########################################

dat <- melt(cls_mx)
colnames(dat) <- c("class","chr","num")

library(ggplot2)

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
            "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
            "#E34234", "#9A2A2A", "#630330", "#342D7E",
            "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

theme_set(
  theme_void() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)


p <- ggplot(dat, aes(x="", y=num, fill=class, group=class))
p + geom_bar(width = 1, stat="identity")+
  coord_polar("y")+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  ggtitle("TE Distribuiton in Human Genome by Class")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "../figure/pie_genome_te.cls_dist.png")


#########################################
###### class bar chart for diff chrs ####
#########################################

### fill

dat$chr = factor(dat$chr, levels = paste0("chr", c(seq(1,22),"X", "Y", "M", "Un")))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)

p <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=class))
p + geom_col(position="fill")+
  labs(title = "TE Class Distribution in Diff. Chromosomes",
       x = NULL,y = 'log10 TE Class Proportion')+
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 8),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()

ggsave(filename = "../figure/genome.class_log_dist.fill.png")
ggsave(filename = "../figure/genome.class_log_dist.fill.pdf")


### stack
#### log10
theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)

q <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=class))
q + geom_col(position="stack")+
  labs(title = "TE Class Distribution in Diff. Chromosomes",
       x = NULL,y = 'log10 TE Class Number')+
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()

ggsave(filename = "../figure/genome.class_log_dist.stack.png")
ggsave(filename = "../figure/genome.class_log_dist.stack.pdf")


### stack
#### no log tranformation
theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)

q <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=class))
q + geom_col(position="stack")+
  labs(title = "TE Class Distribution in Diff. Chromosomes",
       x = NULL,y = 'TE Class Number')+
  coord_flip()+
  scale_fill_manual(values = mycols)+
  guides(fill = guide_legend(title="TE Class"))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave(filename = "../figure/genome.class_dist.stack.png")
ggsave(filename = "../figure/genome.class_dist.stack.pdf")


##########################################
######## construct matrix | family #######
##########################################
fml_mx <- matrix(data=NA, nrow = length(all.fml), ncol = length(genome_te$chr_inc %>% unique()),
                 byrow = FALSE, dimnames = list(all.fml,genome_te$chr_inc %>% unique()))

for (fml in all.fml) {
  for (j in 1:length(genome_te$chr_inc %>% unique())) {
    
    fq <- as.data.frame(table(filter(genome_te,
                                     chr_inc==unique(genome_te$chr_inc)[j])[,"family"]))
    
    if (fq %>% nrow() < length(all.fml)) {
      absc <- matrix(nrow = length(setdiff(all.fml,fq$Var1)),
                     ncol = 2)
      absc[,1] <-  setdiff(all.fml,fq$Var1)
      absc[,2] <- rep(0, length(setdiff(all.fml,fq$Var1)))
      colnames(absc) <- c("Var1", "Freq")
      fq <- rbind(fq,absc)
    }
    fml_mx[fml,j] <- (fq %>% filter(Var1 == fml))[,2]
  }
}


#################=============== plot ================#########################


#############################################
######## family pie chart for genome ########
#############################################
dat <- melt(fml_mx)
colnames(dat) <- c("fml","chr","num")

library(ggplot2)


theme_set(
  theme_void() +
    theme(legend.position = "bottom",
          legend.key.size = unit(12, "pt"))
)

mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
             "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
             "#E34234", "#9A2A2A", "#630330", "#342D7E",
             "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

morecols <- (grDevices::colorRampPalette(mycols))(56) ### expand the platte!

p <- ggplot(dat, aes(x="", y=num, fill=fml, group=fml))
p + geom_bar(width = 1, stat="identity")+
  coord_polar("y")+
  scale_fill_manual(values = morecols)+ 
  guides(fill = guide_legend(title="TE Family"))+
  ggtitle("TE Distribuiton in Human Genome by Family")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.margin=unit(c(1.5, 1, 1.5, 1),'lines'))
ggsave(filename = "../figure/pie_genome_te.fml_dist.png", width = 7, height = 6)


#########################################
##### family bar chart for diff chrs ####
#########################################

### fill

dat$chr = factor(dat$chr, levels = paste0("chr", c(seq(1,22),"X", "Y", "M", "Un")))

theme_set(
  theme_minimal() +
    theme(legend.position = "right",
          legend.key.size = unit(12, "pt"))
)


mycols <- c("#D6EAF8", "#A9DFBF", "#74c493", "#16A085",
            "#aacc22", "#FFBF00", "#E4BF80", "#DC7633",
            "#E34234", "#9A2A2A", "#630330", "#342D7E",
            "#4863A0", "#98AFC7", "#BCC6CC", "#CFD8DC") 

morecols <- (grDevices::colorRampPalette(mycols))(56) ### expand the platte!


p <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=fml))
p + geom_col(position="fill")+
  labs(title = "TE Family Distribution in Diff. Chromosomes",
       x = NULL,y = 'log10 TE Family Proportion')+
  coord_flip()+
  scale_fill_manual(values = morecols)+
  guides(fill = guide_legend(title="TE Family",ncol = 2))+
  theme(axis.text.x = NULL,
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()

ggsave(filename = "../figure/genome.family_log_dist.fill.png", width = 7, height = 6)
ggsave(filename = "../figure/genome.family_log_dist.fill.pdf", width = 7, height = 6)


### stack
q <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=fml))
q + geom_col(position="stack")+
  labs(title = "TE Family Distribution in Diff. Chromosomes",
       x = NULL,y = 'log10 TE Family Number')+
  coord_flip()+
  scale_fill_manual(values = morecols)+
  guides(fill = guide_legend(title="TE Family",ncol = 2))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()

ggsave(filename = "../figure/genome.family_log_dist.stack.png", width = 7, height = 6)
ggsave(filename = "../figure/genome.family_log_dist.stack.pdf", width = 7, height = 6)


### stack
#### no log
q <- ggplot(data = dat %>% filter(chr != "chrM")
            , aes(x=chr, y=as.numeric(num), fill=fml))
q + geom_col(position="stack")+
  labs(title = "TE Family Distribution in Diff. Chromosomes",
       x = NULL,y = 'TE Family Number')+
  coord_flip()+
  scale_fill_manual(values = morecols)+
  guides(fill = guide_legend(title="TE Family",ncol = 2))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))+
  theme(plot.title = element_text(hjust = 0.5))


ggsave(filename = "../figure/genome.family_dist.stack.png", width = 7, height = 6)
ggsave(filename = "../figure/genome.family_dist.stack.pdf", width = 7, height = 6)

