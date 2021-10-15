TWASEWAS <- function(sig_gene_list = NA, ewas_species = NA, twas_species = NA, twas_category = NA, 
                     twas_database = NA, ewas_study = NA, topX = NA, directory = NA, annotation_file_name = NA, num_permutation = 100,
                     bg_file_dir = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/bg_AllSpecies_Amin.RData"){
  
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func3_get_perm_pvalue.R")
  setwd(directory)

  n_top = nrow(pos)
  
  ##############
  # Check database and annotation file are consistent
  ##############
  check_database = unlist(str_split(directory, "/"))[grepl("DB_", unlist(str_split(directory, "/")))] %>% gsub(".*_", "", .)
  check_annotation = annotation_file_name %>% gsub("DataAnnotation.csv", "", .)
  if(check_database != check_annotation){
    stop("Error: Annotation file should be consistent with the database")
  }
  
  ##############
  # Load TWAS annotation file
  ##############
  twas.anno=read.csv(annotation_file_name)
  twas.anno <- twas.anno %>% dplyr::mutate(data = paste("Reference_Database/", Index, "_", Reference, ".csv.gz", sep = ""))
  twas.anno <- twas.anno %>% dplyr::filter(Organism != "Nothobranchius_furzeri") # For now remove killifish

  ##############
  # Apply twas filters
  ##############
  # twas_species <- c("Homo_sapiens", "Rattus_norvegicus", "ALl")
  if(is.na(twas_species)||toupper("all") %in% toupper(twas_species)){ # species (Homo_sapiens, Rattus_norvegicus, Mus_musculus, Macaca_fascicularis, all)
    twas.anno = twas.anno
  }else{
    twas.anno <- twas.anno %>% dplyr::filter(toupper(Organism) %in% toupper(twas_species))
  }
  # twas_category = "Intervention"
  if(is.na(twas_category)||toupper("all") %in% toupper(twas_category)){ # category (Age, Intervention, Transcription)
    twas.anno = twas.anno
  }else{
    twas.anno <- twas.anno %>% dplyr::filter(toupper(Category) %in% toupper(twas_category))
  }
  # twas_database = c("Pubmed", "Tabula muris Senis")
  if(is.na(twas_database)||toupper("all") %in% toupper(twas_database)){ # database (Pubmed, GenAge, GEO, Tabula Muris Senis)
    twas.anno = twas.anno
  }else{
    twas.anno <- twas.anno %>% dplyr::filter(toupper(Database) %in% toupper(twas_database))
  }
  print(paste("Number of gene sets used from the database:", nrow(twas.anno)))
  
  ##############
  # Load background, based on Amin's annotation file
  ##############
  load(bg_file_dir)
  print(paste("Now background species include:", paste(names(bg_AllSpecies_Amin), collapse = ", ")))
  # Load ortholog map
  # ortho <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/ncbi_ortho_human_wide_final.csv", header = T)
  # map <- bg %>%
  #   dplyr::left_join(ortho[,c("Homo_sapiens_symbol", "Mus_musculus_symbol")], by = c("SYMBOL" = "Homo_sapiens_symbol"))
  
  ##############
  # Enrichment analysis
  ##############
  print("Start running enrichment analysis")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = nrow(twas.anno), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  output.all={}
  for(k in 1:nrow(twas.anno)){
    tryCatch({
      # print(twas.anno$Reference[k])
      other=read.csv(gzfile(twas.anno$data[k])) # each gene list
      other.name=twas.anno$Trait[k] # trait of this list
      other.index = twas.anno$Index[k] # index of this list
      
      if(twas.anno$Organism[k] == "Homo_sapiens"){ # need to decide which bg should be used based on species, and also restrict the input list
        bg = bg_AllSpecies_Amin[["human"]]
      }else if(twas.anno$Organism[k] == "Rattus_norvegicus"){
        bg = bg_AllSpecies_Amin[["rat"]]
      }else if(twas.anno$Organism[k] == "Mus_musculus"){
        bg = bg_AllSpecies_Amin[["mouse"]]
      }else if(twas.anno$Organism[k] == "Macaca_fascicularis"){
        bg = bg_AllSpecies_Amin[["macaque"]]
      }else{
        bg = bg_AllSpecies_Amin[["human"]]
      }
      
      if(ewas_species == "mouse"){ # if the ewas was based on mouse or rat, need to further restrict the background to the overlap between mouse/rat and TWAS species; need to restrict the gene set as well but will need a ortholog map
        bg <- bg %>% dplyr::filter(CGid %in% bg_AllSpecies_Amin[["mouse"]]$CGid)
      }else if(ewas_species == "rat"){
        bg <- bg %>% dplyr::filter(CGid %in% bg_AllSpecies_Amin[["rat"]]$CGid)
      }
      
      sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg$CGid)
      other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL))) # restrict pathway based on background
      
      # conduct hypergeometric test
      temp_result <- hypercalc(background = bg, target = sig_gene_list_used, pathway = other, Index = other.index)
      output.all <- rbind(output.all, temp_result)
      
      setTxtProgressBar(pb, k)
    },
    error=function(err) {
      print(paste("Error at", twas.anno$Reference[k], "- Skip it"))
    }
    )
  }
  close(pb) # Close the connection
  print("Done!")
  
  output.all <- 
    output.all %>%
    as.data.frame()
  output.all$Index <- as.numeric(output.all$Index)
  output.all <- 
    output.all %>%
    dplyr::left_join(twas.anno, by = "Index") %>%
    dplyr::arrange(as.numeric(P_value))
  
  # control for multiple comparison
  output.all$fdr <- p.adjust(output.all$P_value, method = "BH")
  
  output.all$Ewas.Study = ewas_study
  output.all$Top.N.CpG = topX
  
  ##############
  # Permutation test p values
  ##############
  # Calculate a p-value based on permutation test
  # Repeat many times to populate a list of scores. Using maximum likelihood estimation, these scores are modeled as a Gamma distribution (this is the null distribution), and a cumulative distribution function (CDF) is calculated.
  ## https://stackoverflow.com/questions/45536234/how-would-you-fit-a-gamma-distribution-to-a-data-in-r, https://rpubs.com/mpfoley73/459051
  
  output.all <- output.all %>% arrange(Index)
  
  output.all$perm_p_gamma <- get_perm_pvalue(actual_data = output.all, method = "gamma", ewas_species = ewas_species, 
                                             n_top = n_top, n_perm = num_permutation, background_array = bg_AllSpecies_Amin,
                                             anno_input = twas.anno, directory = directory)
  output.all$perm_p_nonpar <- get_perm_pvalue(actual_data = output.all, method = "nonpar", ewas_species = ewas_species, 
                                              n_top = n_top, n_perm = num_permutation, background_array = bg_AllSpecies_Amin,
                                              anno_input = twas.anno, directory = directory)
  
  output.all <- output.all %>% dplyr::arrange(as.numeric(perm_p_gamma))
  
  # Reorder
  output.all <- output.all %>%
    dplyr::select(Ewas.Study, Top.N.CpG, Index, Reference, Hit, list_size, Actual_pct, Exp_pct, P_value, fdr, perm_p_gamma, perm_p_nonpar,
                  Hit_genes, Hit_cpgs, Technology, Organism, Tissue, Cell, Type, Category, Trait, Note, PMID, Database)
  return(output.all)
}