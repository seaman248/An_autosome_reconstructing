atr_genes <-  cGeneTable %>%
  select(atr_ID, atr_scf, atr_start, atr_end, atr_strand)

scf_order_data <- read.csv2('./atr_order.csv')

# Определяем порядок хромосом
chr_order <- list('2L', '2R', '3L', '3R')
#
scf_order <- list()
# Делаем лист из data.frame'ов для каждой хромосомы 
scf_order <- lapply(chr_order, function(i){
  scf_order$i<-as.data.frame(scf_order_data[scf_order_data$chr==i,])
})
# Изначальная таблица больше не нужна
rm(scf_order_data)

# 
