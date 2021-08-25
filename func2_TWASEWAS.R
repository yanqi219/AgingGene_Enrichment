TWASEWAS <- function(sig_gene_list = pos){
  output.all={}
  for(k in 1:nrow(twas.anno)){
    print(twas.anno$Reference[k])
    other=read.csv(gzfile(twas.anno$data[k])) # each gene list
    other.name=twas.anno$Trait[k] # trait of this list
    other.index = twas.anno$Index[k] # index of this list
    
    if(twas.anno$Organism[k] == "Homo_sapiens"){ # need to decide which bg should be used based on species, and also restrict the input list
      bg = bg_human
      sig_gene_list_used = sig_gene_list
    }else if(twas.anno$Organism[k] == "Rattus_norvegicus"){
      bg = bg_rat
      sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_rat$CGid)
    }else if(twas.anno$Organism[k] == "Mus_musculus"){
      bg = bg_mouse
      sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_mouse$CGid)
    }else if(twas.anno$Organism[k] == "Macaca_fascicularis"){
      bg = bg_macaque
      sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_macaque$CGid)
    }else{
      next
    }
    
    # conduct hypergeometric test
    temp_result <- hypercalc(background = bg, target = sig_gene_list_used, pathway = other, Index = other.index)
    output.all <- rbind(output.all, temp_result)
  }
  output.all <- 
    output.all %>%
    as.data.frame() %>%
    dplyr::left_join(twas.anno, by = "Index") %>%
    dplyr::arrange(P_value)
  
  # control for multiple comparison
  output.all$fdr <- p.adjust(output.all$P_value, method = "BH")
  return(output.all)
}