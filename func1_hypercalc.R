hypercalc <- function(background = bg, target = sig_gene_list_used, pathway = other, Index = other.index){
  # http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
  hit_gene <-
    target %>%
    dplyr::pull(SYMBOL) %>%
    unique() %>%
    toupper() %>%
    intersect(toupper(pathway$Gene.Symbol))
  hit <- length(hit_gene) # success-in-sample
  hit = hit-1
  pathway_size = length(unique(pathway$Gene.Symbol)) # success-in-bkgd
  reference = length(unique(background$SYMBOL)) # bkgd
  input = length(unique(target$SYMBOL)) # sample-size
  
  total_pct <- pathway_size/reference*100
  sig_pct <- (hit+1)/input*100
  
  pvalue <- phyper(hit, pathway_size, reference - pathway_size, input, lower.tail=F) 
  
  temp_result <- c(Index, hit+1, pathway_size, sig_pct, total_pct, pvalue, paste(hit_gene, collapse = "; "))
  names(temp_result) <- c("Index", "Hit", "list_size", "Actual_pct", "Exp_pct", "P_value", "Hit_genes")
  return(temp_result)
}
