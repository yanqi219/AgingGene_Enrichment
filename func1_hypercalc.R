hypercalc <- function(background = bg, target = sig_gene_list_used, pathway = other, Index = other.index){
  hit <-
    target %>%
    dplyr::pull(SYMBOL) %>%
    unique() %>%
    toupper() %>%
    intersect(toupper(pathway$Gene.Symbol)) %>%
    length()
  hit = hit-1
  pathway_size = nrow(pathway)
  reference = nrow(background)
  input = nrow(target)
  
  total_pct <- pathway_size/reference*100
  sig_pct <- (hit+1)/input*100
  
  pvalue <- phyper(hit, pathway_size, reference - pathway_size, input, lower.tail=F) 
  
  temp_result <- c(Index, hit+1, pathway_size, sig_pct, total_pct, pvalue)
  names(temp_result) <- c("Index", "Hit", "list_size", "Actual_pct", "Exp_pct", "P_value")
  return(temp_result)
}