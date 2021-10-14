rm(list=ls())
options(StringsAsFactor=F)
library(tidyverse)

setwd('/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/DB_TWAS/')
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/funcX_perm_utils.R")
top_num_of_probe = 500
ewas_species_name = "human"
num_of_perm = 1000

# Load original permutation data
load(file = paste("Permutation/TWASEWAS_Permutation_", ewas_species_name, "_top", top_num_of_probe, ".RData", sep = ""))

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
  dplyr::filter(!(is.na(seqnames)))
bg_rat <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Rattus_norvegicus.rnor_6.0.101.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(seqnames)))
bg_macaque <- read.csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Macaca_fascicularis.macaca_fascicularis_5.0.100.HorvathMammalMethylChip40.v1.csv.gz"), header = T) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(seqnames)))

# Check whether need to update the permutation data
unperm <- check_perm(permutation_dir = "Permutation/",
                     anno = twas.anno, species_name = ewas_species_name, topN = top_num_of_probe)

# Conduct permutation tests for newly added gene sets
perm_result <- do_perm(Target_species = ewas_species_name, topN = top_num_of_probe, nperm = num_of_perm, anno_input = unperm)

new_pvalue <- perm_result[[1]]
new_actpct <- perm_result[[2]]

# Update permutation database
perm_pvalue <- rbind(perm_pvalue, new_pvalue)
perm_actpct <- rbind(perm_actpct, new_actpct)

save(perm_pvalue, perm_actpct, file = paste("Permutation/TWASEWAS_Permutation_", ewas_species_name, "_top", top_num_of_probe, ".RData", sep = ""))
