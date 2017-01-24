library('dplyr')
library('genoPlotR')

dGeneTable <- read.csv2('./data/2chr_alb_atr_sin_gam_dirty.csv', na.strings = c(''), stringsAsFactors = FALSE)
names(dGeneTable) <- c(
  'alb_ID', 'alb_chr', 'alb_start', 'alb_end', 'alb_strand',
  'atr_ID', 'atr_scf', 'atr_start', 'atr_end',
  'fun_ID', 'fun_scf', 'fun_start', 'fun_end',
  'gam_ID', 'gam_chr', 'gam_start', 'gam_end',
  'sin_ID', 'sin_scf', 'sin_start', 'sin_end',
  'steph_ID', 'steph_scf', 'steph_start', 'steph_end'
)

# select only alb, atr and gam columns
dGeneTable <- dGeneTable %>%
  select(alb_ID, alb_chr, alb_start, alb_end, alb_strand, atr_ID, atr_scf, atr_start, atr_end, gam_ID, gam_chr, gam_start, gam_end)
# remove strings with empty cells
dGeneTable <-  na.omit(dGeneTable)

# remove all duplicates
dGeneTable <- dGeneTable[!(duplicated(dGeneTable$alb_ID, fromLast = TRUE) | duplicated(dGeneTable$alb_ID)), ]
dGeneTable <- dGeneTable[!(duplicated(dGeneTable$atr_ID, fromLast = TRUE) | duplicated(dGeneTable$atr_ID)),]
dGeneTable <- dGeneTable[!(duplicated(dGeneTable$gam_ID, fromLast = TRUE) | duplicated(dGeneTable$gam_ID)),]
# remove unkn chromosom of gam
dGeneTable$alb_chr <- as.factor(dGeneTable$alb_chr)
dGeneTable$atr_scf <- as.factor(dGeneTable$atr_scf)
dGeneTable$gam_chr <- as.factor(dGeneTable$gam_chr)

dGeneTable <- dGeneTable[which(dGeneTable$gam_chr !='UNKN'),]

# Change rank
row.names(dGeneTable) <- 1:nrow(dGeneTable)

# write table to get strand from vectorbase biomart
write.csv2(dGeneTable, './data/GeneTableForGetStrand.csv') # using this table for extract strands from biomart

# read data from biomart
atr_strands <- read.csv('./data/atr_strands.csv')
gam_strands <- read.csv('./data/gam_strands.csv')



# paste this data in exsisting order
dGeneTable$atr_strand <- atr_strands$Strand[match(dGeneTable$atr_ID, atr_strands$Gene.stable.ID)]
dGeneTable$gam_strand <- gam_strands$Strand[match(dGeneTable$gam_ID, gam_strands$Gene.stable.ID)]

# reorder columns
dGeneTable <- dGeneTable[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12, 13, 15)]
rm(atr_strands, gam_strands)

# write clean data
#write.csv2(dGeneTable, './data/alb_atr_gam-clean.csv')

#rm(dGeneTable)
