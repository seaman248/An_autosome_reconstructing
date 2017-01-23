atr_genes <-  cGeneTable %>%
  select(atr_ID, atr_scf, atr_start, atr_end, atr_strand)
atr_genes$tStart <- atr_genes$tEnd <- atr_genes$tStrand <- atr_genes$chr <- rep(NA, nrow(atr_genes))
scf_order <- read.csv2('./atr_order.csv')

for(i in 1:nrow(atr_genes)){
  
}