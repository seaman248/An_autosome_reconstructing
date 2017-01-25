blocks <- read.table('./data/blocks.txt')
names(blocks) <- c(
  'ID',
  'alb_chr',
  'alb_start',
  'alb_end',
  'alb_strand',
  
  'gam_chr',
  'gam_start',
  'gam_end',
  'gam_strand',
  
  'atr_chr',
  'atr_start',
  'atr_end',
  'atr_strand'
)

blocks[,c(4, 8, 12)] <- blocks[,c(3, 7, 11)] + blocks[,c(4, 8, 12)]

# generate segs
blocks[blocks$alb_chr == '3R', c(3, 4)] <- 
  blocks[blocks$alb_chr == '3R', c(3, 4)] + max(blocks[blocks$alb_chr == '3L', 4])

b_alb_seq <- dna_seg(data.frame(
  name=as.character(blocks$ID),
  start=blocks$alb_start,
  end=blocks$alb_end,
  strand=blocks$alb_strand,
  gene_type=rep('blocks', nrow(blocks))
))

blocks[blocks$gam_chr == '3L', c(7, 8)] <-
  blocks[blocks$gam_chr == '3L', c(7, 8)] + max(blocks[blocks$gam_chr =='2L', 8])

b_gam_seq <- dna_seg((data_frame(
  name=as.character(blocks$ID),
  start=blocks$gam_start,
  end=blocks$gam_end,
  starnd=blocks$gam_strand,
  gene_types=rep('blocks', nrow(blocks))
)))

b_dna_segs <- list(b_alb_seq)
names(b_dna_segs) <- c('An. albimanus')

plot_gene_map(
  dna_segs = b_dna_segs
)
