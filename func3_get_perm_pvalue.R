get_perm_pvalue <- function(actual_data = output_pos, n_perm = 100, method = "gamma", ewas_species = "human", 
                            n_top, background_array, directory, anno_input){
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_perm.R")
  # Load permutation test results
  if(ewas_species == "human"){
    if(n_top == 50){ # Need to be fixed, 500!!
      load(file = "Permutation/TWASEWAS_Permutation_human_top50.RData")
    }else if(n_top == 500){
      load(file = "Permutation/TWASEWAS_Permutation_human_top500.RData")
    }else if(n_top == 1000){
      load(file = "Permutation/TWASEWAS_Permutation_human_top1000.RData")
    }
  }else if(ewas_species == "mouse"){
    if(n_top == 50){ # Need to be fixed, 500!!
      load(file = "Permutation/TWASEWAS_Permutation_mouse_top50.RData")
    }else if(n_top == 500){
      load(file = "Permutation/TWASEWAS_Permutation_mouse_top500.RData")
    }else if(n_top == 1000){
      load(file = "Permutation/TWASEWAS_Permutation_mouse_top1000.RData")
    }
  }
  
  if(exists("perm_pvalue") & exists("perm_pvalue")){
    print(paste("Have pre calculated permutation data (n_perm = 1000): ", "TWASEWAS_Permutation_", ewas_species, "_top", n_top, sep = ""))
    n_perm = 1000
  }else{
    print("Need to calculate permutation (default n_perm = 100)")
    perm_result <- do_perm(Target_species = ewas_species, topN = n_top, nperm = n_perm, 
                           anno_input = anno_input, background_array = background_array, directory = directory)
    perm_pvalue <- perm_result[[1]]
    perm_actpct <- perm_result[[2]]
  }
  
  perm_pvalue <- perm_pvalue %>% dplyr::filter(Index %in% actual_data$Index) # restrict to actual data
  perm_actpct <- perm_actpct %>% dplyr::filter(Index %in% actual_data$Index)
  
  if(method == "gamma"){
    perm_p <- {}
    for (i in 1:nrow(actual_data)) {
      # temp <- sum(perm_pvalue[i,1:n_perm] <= actual_data$P_value[i])/50n_perm0
      # perm_p <- rbind(perm_p, temp)
      a <- t(perm_pvalue[i,1:n_perm])
      if(round(sum(a),5) < n_perm){              # Some pathways have all 1 p values for permutation tests
        # fitdistrplus::plotdist(as.vector(a), histo = TRUE, demp = TRUE)
        invisible(capture.output(fit.gamma <- fitdistrplus::fitdist(as.vector(a), distr = "gamma", method = "mle")))
        fit.gamma <- summary(fit.gamma)
        temp <- pgamma(q = round(as.numeric(actual_data$P_value[i]),8), shape = fit.gamma$estimate[1], scale = 1/fit.gamma$estimate[2]) # Here I will round the p-value [CHECK]
      }else{
        temp = 1
      }
      perm_p <- rbind(perm_p, temp)
    }
  }else if(method == "nonpar"){
    perm_p <- {}
    for (i in 1:nrow(actual_data)) {
      # temp <- sum(perm_pvalue[i,1:n_perm] <= actual_data$P_value[i])/50n_perm0
      # perm_p <- rbind(perm_p, temp)
      a <- t(perm_pvalue[i,1:n_perm])
      temp <- 1-sum(a > round(as.numeric(actual_data$P_value[i]),8))/n_perm
      perm_p <- rbind(perm_p, temp)
    }
  }
  return(perm_p)
}

