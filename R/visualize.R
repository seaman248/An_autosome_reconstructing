source('./R/2nd clean.R')

# add through coordinate for every sp==============================================================
  #alb---------------------------------------------------------------------------------------------
cGeneTable$alb_start[cGeneTable$alb_chr=='3R'] <- 
  cGeneTable$alb_start[cGeneTable$alb_chr=='3R'] + max(cGeneTable$alb_end[cGeneTable$alb_chr=='3L'])
cGeneTable$alb_end[cGeneTable$alb_chr=='3R'] <- 
  cGeneTable$alb_end[cGeneTable$alb_chr=='3R'] + max(cGeneTable$alb_end[cGeneTable$alb_chr=='3L'])
alb_lims <- c(
  max(cGeneTable$alb_start[cGeneTable$alb_chr=='3L']), 
  min(cGeneTable$alb_end[cGeneTable$alb_chr=='3L']), 
  max(cGeneTable$alb_end[cGeneTable$alb_chr=='3R']),
  min(cGeneTable$alb_start[cGeneTable$alb_chr=='3R'])
  )
  #atr---------------------------------------------------------------------------------------------

atr_lims <- c(
  max(cGeneTable$atr_start[cGeneTable$atr_chr == '2L']),
  min(cGeneTable$atr_end[cGeneTable$atr_chr == '2L']),
  max(cGeneTable$atr_start[cGeneTable$atr_chr == '2R']),
  min(cGeneTable$atr_end[cGeneTable$atr_chr == '2R']),
  min(cGeneTable$atr_start[cGeneTable$atr_chr == '3L']),
  max(cGeneTable$atr_end[cGeneTable$atr_chr == '3L'])
)

  #gam---------------------------------------------------------------------------------------------

cGeneTable$gam_start[cGeneTable$gam_chr=='2L'] <- 
  cGeneTable$gam_start[cGeneTable$gam_chr=='2L'] + max(cGeneTable$gam_end[cGeneTable$gam_chr=='3L'])

cGeneTable$gam_end[cGeneTable$gam_chr=='2L'] <- 
  cGeneTable$gam_end[cGeneTable$gam_chr=='2L'] + max(cGeneTable$gam_end[cGeneTable$gam_chr=='3L'])

gam_lims <- c(
  max(cGeneTable$gam_start[cGeneTable$gam_chr=='2L']),
  min(cGeneTable$gam_end[cGeneTable$gam_chr=='2L']),
  min(cGeneTable$gam_start[cGeneTable$gam_chr=='3L']),
  max(cGeneTable$gam_end[cGeneTable$gam_chr=='3L'])
)

#create comparisons

alb_atr_comparison <- as.comparison(data.frame(
  start1=cGeneTable$alb_start, end1=cGeneTable$alb_end,
  start2=cGeneTable$atr_start, end2=cGeneTable$atr_end
))

atr_gam_comparition <- as.comparison(data.frame(
  start1=cGeneTable$atr_start, end1=cGeneTable$atr_end,
  start2=cGeneTable$gam_start, end2=cGeneTable$gam_end
))

# create dna_seqs
alb_col <- rep('red', nrow(cGeneTable))
alb_col[cGeneTable$alb_strand=='-1'] <- 'blue'

alb_seq <- dna_seg(data.frame(
  name=cGeneTable$alb_ID,
  start=cGeneTable$alb_start,
  end=cGeneTable$alb_end,
  strand=cGeneTable$alb_strand,
  col=alb_col,
  gene_type=rep('blocks', nrow(cGeneTable))
))

atr_col <- rep('red', nrow(cGeneTable))
atr_col[cGeneTable$atr_strand=='-1'] <- 'blue'

atr_seq <- dna_seg(data.frame(
  name=cGeneTable$atr_ID,
  start=cGeneTable$atr_start,
  end=cGeneTable$atr_end,
  strand=as.integer(cGeneTable$atr_strand),
  col=atr_col,
  gene_type=rep('blocks', nrow(cGeneTable))
))

gam_col <- rep('red', nrow(cGeneTable))
gam_col[cGeneTable$gam_strand=='-1'] <- 'blue'

gam_seq <- dna_seg(data.frame(
  name=cGeneTable$gam_ID,
  start=cGeneTable$gam_start,
  end=cGeneTable$gam_end,
  strand=cGeneTable$gam_strand,
  col=gam_col,
  gene_type=rep('blocks', nrow(cGeneTable))
))

dna_seqs <- list(alb_seq, atr_seq, gam_seq)
names(dna_seqs) <- c('alb 3L-\\-3R(rev)', 'atr 3L-\\-2R-\\-2L', 'gam 3L-\\-2L')
#plot
plot_gene_map(
  dna_segs = dna_seqs,
  comparisons = list(alb_atr_comparison, atr_gam_comparition),
  xlims = list(alb_lims, atr_lims, gam_lims)
)
