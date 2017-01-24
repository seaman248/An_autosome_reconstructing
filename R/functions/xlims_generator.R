xlimsGen <- function (chr_order, strand_order, start_v, end_v, chr_v){
  xlims <- lapply(chr_order, function(chr_i, strand, start_v, end_v, chr_v){
    if(strand){
      c(min(start_v[chr_v == chr_i]),
        max(end_v[chr_v == chr_i]))
    } else {
      c(max(end_v[chr_v == chr_i]),
        min(start_v[chr_v == chr_i]))
    }
  }, strand=strand_order, start_v=start_v, end_v=end_v, chr_v=chr_v)
  return(unlist(xlims))
}