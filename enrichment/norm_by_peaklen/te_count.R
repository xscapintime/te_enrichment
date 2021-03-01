# never mind the chromosomes
## for TE count

rm(list = ls())
library(dplyr)
library(reshape2)

options(scipen = 999)

## load rmsk table
rmsk <- read.table("../../repeatmasker_track/del_qmarkline_rmsk.bed", header = F)
genome_te <- cbind(rmsk$V1, rmsk$V2, rmsk$V3, rmsk$V4 , paste(rmsk$V4,rmsk$V5)) %>% as.data.frame() %>% tibble()

### TE class count
te_cls_count <- genome_te %>% 
  count(V4, sort = TRUE) %>% 
  mutate(prop = n / sum(n))
save(te_cls_count,file = "te_cls_count.Rdata")

### TE famliy count
te_fml_count <- genome_te %>% 
  count(V4,V5, sort = TRUE) %>% 
  mutate(prop = n / sum(n))
save(te_fml_count,file = "te_fml_count.Rdata")



