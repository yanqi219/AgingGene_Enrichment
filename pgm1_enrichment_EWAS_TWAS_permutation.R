rm(list=ls())
options(StringsAsFactor=F)
library(tidyverse)

 # Change to different database
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_perm.R")
top_num_of_probe = 1000
ewas_species_name = "human"
num_of_perm = 1000
working_dir <- '/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/DB_TWAS/'
setwd(working_dir)

# Load original permutation data
load(file = paste("Permutation/TWASEWAS_Permutation_", ewas_species_name, "_top", top_num_of_probe, ".RData", sep = ""))

# Load TWAS annotation file
twas.anno=read.csv('TWASDataAnnotation.csv') # Change to different annotation files
twas.anno <- twas.anno %>% dplyr::mutate(data = paste("Reference_Database/", Index, "_", Reference, ".csv.gz", sep = ""))
twas.anno <- twas.anno %>% dplyr::filter(Organism != "Nothobranchius_furzeri") # For now remove killifish

# Load background, based on Amin's annotation file
load("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/bg_AllSpecies_Amin.RData")

# Check whether need to update the permutation data
unperm <- check_perm(permutation_dir = "Permutation/",
                     anno = twas.anno, species_name = ewas_species_name, topN = top_num_of_probe)

# Conduct permutation tests for newly added gene sets
perm_result <- do_perm(Target_species = ewas_species_name, topN = top_num_of_probe, nperm = num_of_perm, 
                       anno_input = unperm, background_array = bg_AllSpecies_Amin, directory = working_dir)

new_pvalue <- perm_result[[1]]
new_actpct <- perm_result[[2]]

# Update permutation database
if(exists("perm_pvalue") & exists("perm_actpct")){
  perm_pvalue <- rbind(perm_pvalue, new_pvalue)
  perm_actpct <- rbind(perm_actpct, new_actpct)
}else{
  perm_pvalue <- new_pvalue
  perm_actpct <- new_actpct
}

save(perm_pvalue, perm_actpct, file = paste("Permutation/TWASEWAS_Permutation_", ewas_species_name, "_top", top_num_of_probe, ".RData", sep = ""))
