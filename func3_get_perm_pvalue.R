get_perm_pvalue <- function(actual_data = output_pos, n_perm = 1000, method = "nonpar", ewas_species = "human"){
  # Load permutation test results
  if(ewas_species == "human"){
    load(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_Permutation_human.RData")
  }else if(ewas_species == "mouse"){
    load(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_Permutation_mouse.RData")
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
        fit.gamma <- fitdistrplus::fitdist(as.vector(a), distr = "gamma", method = "mle")
        fit.gamma <- summary(fit.gamma)
        temp <- pgamma(q = round(as.numeric(actual_data$P_value[i]),2), shape = fit.gamma$estimate[1], scale = 1/fit.gamma$estimate[2]) # Here I will round the p-value [CHECK]
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
      temp <- 1-sum(a > round(as.numeric(actual_data$P_value[i]),2))/n_perm
      perm_p <- rbind(perm_p, temp)
    }
  }
  return(perm_p)
}

