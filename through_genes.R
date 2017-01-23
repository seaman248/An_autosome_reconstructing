atr_genes <-  cGeneTable %>%
  select(atr_ID, atr_scf, atr_start, atr_end, atr_strand)
atr_genes$tStart <- atr_genes$tEnd <- atr_genes$tStrand <- atr_genes$chr <- rep(NA, nrow(atr_genes))
scf_order <- read.csv2('./atr_order.csv')

# Делаем лист из data.frame'ов для каждой хромосомы 
scf_order <- lapply(chr_order, function(i){
  scf_order$i<-scf_order_data[scf_order_data$chr==i,]
})
names(scf_order) <- chr_order
# Изначальная таблица больше не нужна
rm(scf_order_data)

# 


for(i in 1:nrow(atr_genes)){
  
}