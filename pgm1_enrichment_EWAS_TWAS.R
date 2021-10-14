#
# very important: change background
# no chr23 and 24 for other GWAS such as GWAS of Education***very important
#
rm(list=ls())
options(StringsAsFactor=F)
library(tidyverse)

setwd('/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/DB_TWAS/')
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func2_TWASEWAS.R")

################## INPUT PARAMETER ##################
# Load targeted gene list
# Need one column named "SYMBOL" and one column named "CGid"!!

# Max lifespan
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, EWAS max lifspan.csv", header = T) %>%
#   dplyr::filter(group == "positive") %>%
#   dplyr::filter(Tissue == "ALL") # Top 500 CpGs related to max lifespan
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, EWAS max lifspan.csv", header = T) %>%
#   dplyr::filter(group == "negative") %>%
#   dplyr::filter(Tissue == "ALL") # Top 500 CpGs related to max lifespan

# Max lifespan, adj phylo
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, Phylogenetic EWAS of maximum lifespan V2.csv", header = T) %>%
#   dplyr::filter(group == "positive") %>%
#   dplyr::filter(Tissue == "Phylo_ALL") # Top 500 CpGs related to max lifespan
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, Phylogenetic EWAS of maximum lifespan V2.csv", header = T) %>%
#   dplyr::filter(group == "negative") %>%
#   dplyr::filter(Tissue == "Phylo_ALL") # Top 500 CpGs related to max lifespan

# Max lifespan, adj weight
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_allEWAS_weightAdjusted_lifespan.csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta >= 0)
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_allEWAS_weightAdjusted_lifespan.csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta < 0)

# Averaged Maturity
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_phylo_allEWAS_averagedMaturity.yrs.csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta >= 0)
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_phylo_allEWAS_averagedMaturity.yrs.csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta < 0)

# Gestation.Incubation
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_phylo_allEWAS_Gestation.Incubation..days..csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta >= 0)
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/topCGs_EWAS_Caesar/Eutherians_phylo_allEWAS_Gestation.Incubation..days..csv", header = T) %>%
#   dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
#   dplyr::mutate(Meta = as.numeric(Meta)) %>%
#   dplyr::filter(Meta < 0)

# Chronological age
# pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/Top1000PosNegCpGs_EWASofAge_Ake.csv", header = T) %>%
#   dplyr::filter(class == "pos") %>%
#   dplyr::filter(tissue == "All") %>%
#   dplyr::rename(SYMBOL = "Gene") %>%
#   dplyr::rename(CGid = "CpG")
# neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Important_Genes_HorvathLab/Top1000PosNegCpGs_EWASofAge_Ake.csv", header = T) %>%
#   dplyr::filter(class == "neg") %>%
#   dplyr::filter(tissue == "All") %>%
#   dplyr::rename(SYMBOL = "Gene") %>%
#   dplyr::rename(CGid = "CpG")

ewas_species = "human" # currently it can be human, mouse, or rat; will add more later
#####################################################

################## MAIN FUNCTION ##################
output_pos <- TWASEWAS(sig_gene_list = pos, ewas_species = ewas_species, twas_species = NA, twas_category = NA, twas_database = NA) %>%
  dplyr::mutate(Direction = "Hyper")
output_neg <- TWASEWAS(sig_gene_list = neg, ewas_species = ewas_species, twas_species = NA, twas_category = NA, twas_database = NA) %>%
  dplyr::mutate(Direction = "Hypo")
###################################################

output_all <- rbind(output_pos, output_neg) %>%
  arrange(perm_p_nonpar)

write.table(output_all,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_enriched_results_averagedMaturity_phylo.csv",sep=',',row.names = F,quote=F)

################## DRAW PLOT ##################
plot <-                # Plot those with permutation p value < 0.05 and hit > 5
  output_all %>%
  dplyr::filter(perm_p_nonpar < 0.05) %>%
  dplyr::filter(as.numeric(Hit) >= 3) %>%
  dplyr::mutate(p.value.log10 = -log10(perm_p_nonpar+0.001)) %>%
  dplyr::mutate(Actual_pct = as.numeric(Actual_pct)) %>%
  dplyr::mutate(Direction = factor(Direction, levels = c("Hyper", "Hypo")))
order <- unique(plot$Reference)
plot <- plot %>%
  dplyr::mutate(Reference = factor(Reference, levels = rev(order)))

png(file="/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_enriched_results_averagedMaturity_phylo.png", width=1600, height=800)
ggplot(plot, aes_string(x="Actual_pct", y="Reference", size="Actual_pct", color="p.value.log10")) + 
  geom_point() +
  scale_color_continuous(low="blue", high="red", name = "-Log10(P-value)", guide=guide_colorbar(reverse=F)) +
  ylab(NULL) + ggtitle("") + DOSE::theme_dose(24) + scale_size(range=c(6, 16)) +
  facet_grid(. ~ Direction)
dev.off()
################################################