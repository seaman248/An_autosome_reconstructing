source('./R/clean_first_data.R')
source('./R/functions/through_genes.R')

# generate atr through coordinates of genes

atr_genes <-  dGeneTable %>%
  select(atr_ID, atr_scf, atr_start, atr_end, atr_strand)

atr_scf_order <- read.csv2('./data/atr_order.csv')
sum(atr_scf_order$Size)

atr_genes <- through_num(atr_genes, atr_scf_order)
max(na.omit(atr_genes2$tEnd))

# insert atr_genes into GENETABLE

cGeneTable <- dGeneTable[,c(1:5, 11:15)]
rm(dGeneTable)
cGeneTable[,11:15] <- atr_genes
names(cGeneTable) <- c(names(cGeneTable[,1:10]), 'atr_ID', 'atr_chr', 'atr_start', 'atr_end', 'atr_strand')
rm(atr_genes)
# clean data
cGeneTable <- na.omit(cGeneTable)

cGeneTable <- cGeneTable %>%
  filter(gam_chr=='3L'|gam_chr=='2L', atr_chr!='3R')
