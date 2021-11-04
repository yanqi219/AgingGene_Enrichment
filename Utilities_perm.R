####################
# Check whether need to update the permutation data
####################
check_perm <- function(permutation_dir, anno, species_name, topN){
  library(tidyr)
  permutation_file <- list.files(permutation_dir, pattern = ".RData")
  permutation_file <- permutation_file[grepl(paste(species_name, "_top", topN, ".RData", sep = ""), permutation_file)]
  if(length(permutation_file)==0){
    print("Don't have this permutation data")
    unperm <- anno
    print(paste("Need to update", nrow(unperm), "new records"))
  }else{
    load(paste(permutation_dir, permutation_file, sep = ""))
    if(nrow(perm_pvalue) != nrow(anno)){
      unperm <- anno %>% dplyr::filter(!(Index %in% perm_pvalue$Index))
      print(paste(permutation_file, "need to update", nrow(unperm), "new records"))
    }else{
      print(paste(permutation_file, "is up-to-date"))
    }
  }
  return(unperm)
}

####################
# TWASEWAS_perm function
####################
TWASEWAS_perm <- function(sig_gene_list = pos, ewas_species = ewas_species, twas_species = NA, twas_category = NA, twas_database = NA, anno, background_array, directory){
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
  setwd(directory)
  # Apply twas filters
  # twas_species <- c("Homo_sapiens", "Rattus_norvegicus", "ALl")
  if(is.na(twas_species)||toupper("all") %in% toupper(twas_species)){ # species (Homo_sapiens, Rattus_norvegicus, Mus_musculus, Macaca_fascicularis, all)
    anno = anno
  }else{
    anno <- anno %>% dplyr::filter(toupper(Organism) %in% toupper(twas_species))
  }
  # twas_category = "InterVention"
  if(is.na(twas_category)||toupper("all") %in% toupper(twas_category)){ # category (Age, Intervention, Transcription)
    anno = anno
  }else{
    anno <- anno %>% dplyr::filter(toupper(Category) %in% toupper(twas_category))
  }
  # twas_database = c("Pubmed", "Tabula muris Senis")
  if(is.na(twas_database)||toupper("all") %in% toupper(twas_database)){ # database (Pubmed, GenAge, GEO, Tabula Muris Senis)
    anno = anno
  }else{
    anno <- anno %>% dplyr::filter(toupper(Database) %in% toupper(twas_database))
  }
  
  output.all={}
  for(k in 1:nrow(anno)){
    # print(anno$Reference[k])
    other=read.csv(gzfile(anno$data[k])) # each gene list
    other.name=anno$Trait[k] # trait of this list
    other.index = anno$Index[k] # index of this list
    
    if(anno$Organism[k] == "Homo_sapiens"){ # need to decide which bg should be used based on species, and also restrict the input list
      bg = background_array[["human"]]
      # sig_gene_list_used = sig_gene_list
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL))) 
    }else if(anno$Organism[k] == "Rattus_norvegicus"){
      bg = background_array[["rat"]]
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_rat$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
    }else if(anno$Organism[k] == "Mus_musculus"){
      bg = background_array[["mouse"]]
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_mouse$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
    }else if(anno$Organism[k] == "Macaca_fascicularis"){
      bg = background_array[["macaque"]]
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_macaque$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
    }else if(anno$Organism[k] == "Macaca_mulatta"){
      bg = background_array[["rhesus"]]
    }else{
      bg = background_array[["human"]]
    }
    
    if(ewas_species == "mouse"){ # if the ewas was based on mouse or rat, need to further restrict the background to the overlap between mouse/rat and TWAS species; need to restrict the gene set as well but will need a ortholog map
      bg <- bg %>% dplyr::filter(CGid %in% background_array[["mouse"]]$CGid)
    }else if(ewas_species == "rat"){
      bg <- bg %>% dplyr::filter(CGid %in% background_array[["rat"]]$CGid)
    }
    
    sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg$CGid)
    other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL))) # restrict pathway based on background
    
    # conduct hypergeometric test
    temp_result <- hypercalc(background = bg, target = sig_gene_list_used, pathway = other, Index = other.index)
    output.all <- rbind(output.all, temp_result)
  }
  output.all <- 
    output.all %>%
    as.data.frame()
  output.all$Index <- as.numeric(output.all$Index)
  output.all <- 
    output.all %>%
    dplyr::left_join(anno, by = "Index") %>%
    dplyr::arrange(P_value)
  
  # control for multiple comparison
  output.all$fdr <- p.adjust(output.all$P_value, method = "BH")
  return(output.all)
}

####################
# Do permutation test
####################
do_perm <- function(Target_species, topN, nperm, anno_input, background_array, directory){
  setwd(directory)
  start <- Sys.time()
  perm_pvalue = {}
  perm_actpct = {}
  for (i in 1:nperm) {
    print(paste("Permutation round", i))
    set.seed(1106+i)
    perm <- # Randomlly sample 500 rows/1000 rows/50 rows/100 rows
      dplyr::sample_n(background_array[["human"]], topN)
    temp_perm <- TWASEWAS_perm(sig_gene_list = perm, ewas_species = Target_species, anno = anno_input, background_array = background_array, directory = directory) %>% dplyr::arrange(Index)
    temp_pvalue <- as.numeric(temp_perm$P_value)
    temp_actpct <- as.numeric(temp_perm$Actual_pct)
    
    perm_pvalue <- cbind(perm_pvalue, as.numeric(temp_pvalue))
    perm_actpct <- cbind(perm_actpct, as.numeric(temp_actpct))
  }
  perm_pvalue <- perm_pvalue %>% as.data.frame() %>% dplyr::mutate(Reference = anno_input$Reference)
  perm_actpct <- perm_actpct %>% as.data.frame() %>% dplyr::mutate(Reference = anno_input$Reference)
  perm_pvalue$Index <- anno_input$Index
  perm_actpct$Index <- anno_input$Index
  print(Sys.time()-start)
  
  perm_final <- list(perm_pvalue, perm_actpct)
  return(perm_final)
}
