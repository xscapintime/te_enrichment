g.te_cls <- readRDS("../teclass_mx_bychr.rds")
order1 <- row.names(g.te_cls)
g.te_cls <- apply(g.te_cls,2,as.numeric)
row.names(g.te_cls) <- order1

class_order <- apply(g.te_cls,1,sum) %>% sort(decreasing  = T) %>% names()
save(class_order,file = "cls_order_bycontent.Rdata")
