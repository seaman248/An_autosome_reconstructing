

levels(dGeneTable$atr_scf)


# Сквозная нумерация на основе порядка скэфолдов
TRUE_scf_ORDER <- as.factor(sample(levels(dGeneTable$atr_scf)))

atr_genes <- dGeneTable %>% 
  select(atr_ID, atr_scf, atr_start, atr_end, atr_strand)

start_n <- 0

for(i in 1:length(TRUE_scf_ORDER)){
  poss <- which(!is.na(match(atr_genes$atr_scf, TRUE_scf_ORDER[i])))
  atr_genes$atr_start[poss] <- atr_genes$atr_start[poss] + start_n
  atr_genes$atr_end[poss] <- atr_genes$atr_end[poss] + start_n
  start_n <- start_n + max(atr_genes$atr_end[poss])
  rm(poss)
}

