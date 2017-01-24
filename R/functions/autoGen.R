annot_gen <- function (sp_lims, sp_chr_order){
  sp_annot <- annotation(
    x1=sp_lims[seq(1, length(sp_lims), 2)],
    x2=sp_lims[seq(2, length(sp_lims), 2)],
    text=c(sp_chr_order)
  )
}