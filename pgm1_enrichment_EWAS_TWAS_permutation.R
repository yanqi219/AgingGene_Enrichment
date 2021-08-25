rm(list=ls())
options(StringsAsFactor=F)
library(tidyverse)

setwd('/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/')
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func1_hypercalc.R")
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func2_TWASEWAS.R")

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

# Load ortholog map
# ortho <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/ncbi_ortho_human_wide_final.csv", header = T)
# map <- bg %>%
#   dplyr::left_join(ortho[,c("Homo_sapiens_symbol", "Mus_musculus_symbol")], by = c("SYMBOL" = "Homo_sapiens_symbol"))

# Generate random targeted gene list
start <- Sys.time()
perm_pvalue = {}
perm_actpct = {}
for (i in 1:1000) {
  set.seed(1106+i)
  perm <- # Randomlly sample 500 rows
    dplyr::sample_n(bg_human, 500)
  temp_perm <- TWASEWAS(sig_gene_list = perm) %>% dplyr::arrange(Index)
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

save(perm_pvalue, perm_actpct, file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_Permutation.RData")
