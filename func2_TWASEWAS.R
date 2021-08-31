TWASEWAS <- function(sig_gene_list = pos, ewas_species = ewas_species, twas_species = NA, twas_category = NA, twas_database = NA){
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
  
  ##############
  # Load TWAS annotation file
  ##############
  twas.anno=read.csv('TWASDataAnnotation.csv')
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
  # twas_category = "InterVention"
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
  
  ##############
  # Load background, based on Amin's annotation file
  ##############
  bg_human <- readRDS("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Human.Homo_sapiens.hg19.Amin.V7.RDS") %>%
    dplyr::filter(!grepl("rs", CGid)) %>%
    dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid)))
  bg_mouse <- readRDS("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Mouse.Mus_musculus.GRCm38.100.Amin.V6.RDS") %>%
    dplyr::filter(!grepl("rs", CGid)) %>%
    dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
    dplyr::filter(!(is.na(SYMBOL)))
  bg_rat <- readRDS("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Rat.Rattus_norvegicus.Rnor_6.0.101.Amin.V2.RDS") %>%
    dplyr::filter(!grepl("rs", CGid)) %>%
    dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
    dplyr::filter(!(is.na(SYMBOL)))
  bg_macaque <- readRDS("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Macaca_fascicularis.Macaca_fascicularis_5.0.100.Mingjia.AminV2.RDS") %>%
    dplyr::filter(!grepl("rs", CGid)) %>%
    dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
    dplyr::filter(!(is.na(SYMBOL)))
  
  # Load ortholog map
  # ortho <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/ncbi_ortho_human_wide_final.csv", header = T)
  # map <- bg %>%
  #   dplyr::left_join(ortho[,c("Homo_sapiens_symbol", "Mus_musculus_symbol")], by = c("SYMBOL" = "Homo_sapiens_symbol"))
  
  ##############
  # Enrichment analysis
  ##############
  output.all={}
  for(k in 1:nrow(twas.anno)){
    print(twas.anno$Reference[k])
    other=read.csv(gzfile(twas.anno$data[k])) # each gene list
    other.name=twas.anno$Trait[k] # trait of this list
    other.index = twas.anno$Index[k] # index of this list
    
    if(twas.anno$Organism[k] == "Homo_sapiens"){ # need to decide which bg should be used based on species, and also restrict the input list
      bg = bg_human
    }else if(twas.anno$Organism[k] == "Rattus_norvegicus"){
      bg = bg_rat
    }else if(twas.anno$Organism[k] == "Mus_musculus"){
      bg = bg_mouse
    }else if(twas.anno$Organism[k] == "Macaca_fascicularis"){
      bg = bg_macaque
    }else{
      next
    }
    
    if(ewas_species == "mouse"){ # if the ewas was based on mouse or rat, need to further restrict the background to the overlap between mouse/rat and TWAS species; need to restrict the gene set as well but will need a ortholog map
      bg <- bg %>% dplyr::filter(CGid %in% bg_mouse$CGid)
    }else if(ewas_species == "rat"){
      bg <- bg %>% dplyr::filter(CGid %in% bg_rat$CGid)
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
    dplyr::left_join(twas.anno, by = "Index") %>%
    dplyr::arrange(as.numeric(P_value))
  
  # control for multiple comparison
  output.all$fdr <- p.adjust(output.all$P_value, method = "BH")
  
  ##############
  # Permutation test p values
  ##############
  # Calculate a p-value based on permutation test
  # Repeat many times to populate a list of scores. Using maximum likelihood estimation, these scores are modeled as a Gamma distribution (this is the null distribution), and a cumulative distribution function (CDF) is calculated.
  ## https://stackoverflow.com/questions/45536234/how-would-you-fit-a-gamma-distribution-to-a-data-in-r, https://rpubs.com/mpfoley73/459051
  
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func3_get_perm_pvalue.R")
  output.all <- output.all %>% arrange(Index)
  
  output.all$perm_p_gamma <- get_perm_pvalue(actual_data = output.all, method = "gamma", ewas_species = ewas_species)
  output.all$perm_p_nonpar <- get_perm_pvalue(actual_data = output.all, method = "nonpar", ewas_species = ewas_species)
  
  output.all <- output.all %>% dplyr::arrange(as.numeric(perm_p_gamma))
  
  return(output.all)
}