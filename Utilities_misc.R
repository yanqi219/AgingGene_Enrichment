###############################
# Create background R data based on Amin's annotation file - needs to be constantly updated
##############################
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
bg_rhesus <- read_csv(gzfile("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/Mammals/Macaca_mulatta.mmul_10.100.HorvathMammalMethylChip40.v1.csv.gz")) %>%
  dplyr::filter(!grepl("rs", CGid)) %>%
  dplyr::filter(!(is.na(SYMBOL) & !grepl("cg", CGid))) %>%
  dplyr::filter(!(is.na(SYMBOL)))

bg_AllSpecies_Amin <- list(human = bg_human, mouse = bg_mouse, rat = bg_rat, macaque = bg_macaque, rhesus = bg_rhesus)
save(bg_AllSpecies_Amin, file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/bg_AllSpecies_Amin.RData")

###############################
# Reformate Caesar's EWAS output
##############################
library(tidyr)
zscore <- readRDS("/Users/qiyan/Dropbox/Horvath_Lab/MammalianMethCombined/StuffCaesar/JunoFork/allEWAS/Oct2021/Eutherians_allPhyloEWAS_weightAdjusted_zvalue.RDS")
load("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/bg_AllSpecies_Amin.RData")
human_anno <- bg_AllSpecies_Amin[["human"]]
zscore <- zscore %>%
  tibble::rownames_to_column(var = "CGid") %>%
  dplyr::left_join(human_anno, by = "CGid") %>%
  dplyr::arrange(desc(abs(`57.PGLS.EWAS.Log.maxAgeCaesar.OrderALL.TissueALL.N215`))) %>%
  dplyr::mutate(rank_OrderALL_TissueALL = 1:nrow(zscore))

write.csv(zscore, file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/EWASmaxlifespan_topGene_local/topCG_Caesar/Eutherians_allPhyloEWAS_weightAdjusted_zvalue_OCTqy.csv", row.names = F, col.names = T)
