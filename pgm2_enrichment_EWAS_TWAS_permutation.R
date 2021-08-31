rm(list=ls())
options(StringsAsFactor=F)
library(tidyverse)

setwd('/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/')
# source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
# source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func2_TWASEWAS.R")

# Load TWAS annotation file
twas.anno=read.csv('TWASDataAnnotation.csv')
twas.anno <- twas.anno %>% dplyr::mutate(data = paste("Reference_Database/", Index, "_", Reference, ".csv.gz", sep = ""))
twas.anno <- twas.anno %>% dplyr::filter(Organism != "Nothobranchius_furzeri") # For now remove killifish

# Load background, based on Amin's annotation file
bg_human <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Homo_sapiens.hg19.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid)))
bg_mouse <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Mus_musculus.grcm38.100.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(SYMBOL)))
bg_rat <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Rattus_norvegicus.rnor_6.0.101.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(SYMBOL)))
bg_macaque <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Macaca_fascicularis.macaca_fascicularis_5.0.100.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(SYMBOL)))

# Load TWASEWAS_perm function
TWASEWAS_perm <- function(sig_gene_list = pos, ewas_species = ewas_species, twas_species = NA, twas_category = NA, twas_database = NA){
  source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
  
  # Apply twas filters
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
  
  output.all={}
  for(k in 1:nrow(twas.anno)){
    print(twas.anno$Reference[k])
    other=read.csv(gzfile(twas.anno$data[k])) # each gene list
    other.name=twas.anno$Trait[k] # trait of this list
    other.index = twas.anno$Index[k] # index of this list
    
    if(twas.anno$Organism[k] == "Homo_sapiens"){ # need to decide which bg should be used based on species, and also restrict the input list
      bg = bg_human
      # sig_gene_list_used = sig_gene_list
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL))) 
    }else if(twas.anno$Organism[k] == "Rattus_norvegicus"){
      bg = bg_rat
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_rat$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
    }else if(twas.anno$Organism[k] == "Mus_musculus"){
      bg = bg_mouse
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_mouse$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
    }else if(twas.anno$Organism[k] == "Macaca_fascicularis"){
      bg = bg_macaque
      # sig_gene_list_used <- sig_gene_list %>% dplyr::filter(CGid %in% bg_macaque$CGid)
      # other <- other %>% dplyr::filter(toupper(Gene.Symbol) %in% toupper(unique(bg$SYMBOL)))
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
    dplyr::arrange(P_value)
  
  # control for multiple comparison
  output.all$fdr <- p.adjust(output.all$P_value, method = "BH")
  return(output.all)
}

# Generate random targeted gene list for human and mouse
Target_species = "human"

start <- Sys.time()
perm_pvalue = {}
perm_actpct = {}
for (i in 1:1000) {
  set.seed(1106+i)
  perm <- # Randomlly sample 500 rows
    dplyr::sample_n(bg_human, 500)
  temp_perm <- TWASEWAS_perm(sig_gene_list = perm, ewas_species = Target_species) %>% dplyr::arrange(Index)
  temp_pvalue <- as.numeric(temp_perm$P_value)
  temp_actpct <- as.numeric(temp_perm$Actual_pct)
  
  perm_pvalue <- cbind(perm_pvalue, as.numeric(temp_pvalue))
  perm_actpct <- cbind(perm_actpct, as.numeric(temp_actpct))
}
perm_pvalue <- perm_pvalue %>% as.data.frame() %>% dplyr::mutate(Reference = twas.anno$Reference)
perm_actpct <- perm_actpct %>% as.data.frame() %>% dplyr::mutate(Reference = twas.anno$Reference)
perm_pvalue$Index <- twas.anno$Index
perm_actpct$Index <- twas.anno$Index
Sys.time()-start

save(perm_pvalue, perm_actpct, file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_Permutation_human.RData")
