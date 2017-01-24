source('./R/2nd clean.R')
source('./R/functions/xlimsGen.R')
source('./R/functions/autoGen.R')
# add through coordinate for every sp==============================================================
  #alb---------------------------------------------------------------------------------------------
cGeneTable$alb_start[cGeneTable$alb_chr=='3R'] <- 
  cGeneTable$alb_start[cGeneTable$alb_chr=='3R'] + max(cGeneTable$alb_end[cGeneTable$alb_chr=='3L'])
cGeneTable$alb_end[cGeneTable$alb_chr=='3R'] <- 
  cGeneTable$alb_end[cGeneTable$alb_chr=='3R'] + max(cGeneTable$alb_end[cGeneTable$alb_chr=='3L'])

alb_chr_order <- c('3L', '3R')
alb_chr_strand <- c(FALSE, FALSE)

alb_lims <- xlimsGen(alb_chr_order, alb_chr_strand, cGeneTable$alb_start, cGeneTable$alb_end, cGeneTable$alb_chr)
alb_annot <- annot_gen(alb_lims, alb_chr_order)

atr_2R_genes <- cGeneTable[cGeneTable$atr_chr=='2R',]
#alb_annot <- annotation(x1=atr_2R_genes$alb_start, text = atr_2R_genes$alb_ID, rot=30)
  #atr---------------------------------------------------------------------------------------------

atr_chr_order <- c('2L', '2R', '3L')
atr_chr_strand <- c(FALSE, FALSE, TRUE)

atr_lims <- xlimsGen(atr_chr_order, atr_chr_strand, cGeneTable$atr_start, cGeneTable$atr_end, cGeneTable$atr_chr)
atr_annot <- annot_gen(atr_lims, atr_chr_order)

  #gam---------------------------------------------------------------------------------------------

cGeneTable$gam_start[cGeneTable$gam_chr=='2L'] <- 
  cGeneTable$gam_start[cGeneTable$gam_chr=='2L'] + max(cGeneTable$gam_end[cGeneTable$gam_chr=='3L'])

cGeneTable$gam_end[cGeneTable$gam_chr=='2L'] <- 
  cGeneTable$gam_end[cGeneTable$gam_chr=='2L'] + max(cGeneTable$gam_end[cGeneTable$gam_chr=='3L'])

gam_chr_order <- c('2L', '3L')
gam_chr_strand <- c(FALSE, FALSE)

gam_lims <- xlimsGen(gam_chr_order, gam_chr_strand, cGeneTable$gam_start, cGeneTable$gam_end, cGeneTable$gam_chr)
gam_annot <- annot_gen(gam_lims, gam_chr_order)
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
names(dna_seqs) <- c('An. albimanus', 'An. atroparvus', 'An. gambiae')
#plot
plot_gene_map(
  dna_segs = dna_seqs,
  annotations = list(alb_annot, atr_annot, gam_annot),
  comparisons = list(alb_atr_comparison, atr_gam_comparition),
  xlims = list(alb_lims, atr_lims, gam_lims),
  scale = FALSE, dna_seg_scale = TRUE,
  gene_type = 'side_blocks', annotation_height = 2, annotation_cex = 1
)
