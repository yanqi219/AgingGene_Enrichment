#
# very important: change background
# no chr23 and 24 for other GWAS such as GWAS of Education***very important
#
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

# Load targeted gene list
pos <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, EWAS max lifspan.csv", header = T) %>%
    dplyr::filter(group == "positive") %>%
    dplyr::filter(Tissue == "ALL") # Top 500 CpGs related to max lifespan
neg <- read.csv(file = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/Important_genes/Aging_EWAS_HorvathLab/Top hits, EWAS max lifspan.csv", header = T) %>%
  dplyr::filter(group == "negative") %>%
  dplyr::filter(Tissue == "ALL") # Top 500 CpGs related to max lifespan

output_pos <- TWASEWAS(sig_gene_list = pos) %>%
  dplyr::mutate(Direction = "Hyper")
output_neg <- TWASEWAS(sig_gene_list = neg) %>%
  dplyr::mutate(Direction = "Hypo")

# Calculate a p-value based on permutation test
# Repeat many times to populate a list of scores. Using maximum likelihood estimation, these scores are modeled as a Gamma distribution (this is the null distribution), and a cumulative distribution function (CDF) is calculated.
load(file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_Permutation.RData")
n_perm = 1000
output_pos <- output_pos %>% arrange(Index)
output_neg <- output_neg %>% arrange(Index)

perm_p <- {}
for (i in 1:nrow(output_pos)) {
  # temp <- sum(perm_pvalue[i,1:n_perm] <= output_pos$P_value[i])/50n_perm0
  # perm_p <- rbind(perm_p, temp)
  a <- t(perm_pvalue[i,1:n_perm])
  if(sum(a) < n_perm){              # Some pathways have all 1 p values for permutation tests
    fit.gamma <- fitdistrplus::fitdist(as.vector(a), distr = "gamma", method = "mle")
    fit.gamma <- summary(fit.gamma)
    temp <- pgamma(q = output_pos$P_value[i], shape = fit.gamma$estimate[1], scale = 1/fit.gamma$estimate[2]) 
  }else{
    temp = 1
  }
  perm_p <- rbind(perm_p, temp)
}
output_pos$perm_p <- perm_p

perm_p <- {}
for (i in 1:nrow(output_neg)) {
  # temp <- sum(perm_pvalue[i,1:n_perm] <= output_neg$P_value[i])/n_perm
  # perm_p <- rbind(perm_p, temp)
  a <- t(perm_pvalue[i,1:n_perm])
  if(sum(a) < n_perm){              # Some pathways have all 1 p values for permutation tests
    fit.gamma <- fitdistrplus::fitdist(as.vector(a), distr = "gamma", method = "mle")
    fit.gamma <- summary(fit.gamma)
    temp <- pgamma(q = output_neg$P_value[i], shape = fit.gamma$estimate[1], scale = 1/fit.gamma$estimate[2]) 
  }else{
    temp = 1
  }
  perm_p <- rbind(perm_p, temp)
}
output_neg$perm_p <- perm_p

output_all <- rbind(output_pos, output_neg) %>%
  arrange(perm_p)

write.table(output_all,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_enriched_results.csv",sep=',',row.names = F,quote=F)

# Draw plot
plot <-                # Plot those with oermutation p value < 0.05 and hit > 3
  output_all %>%
  dplyr::filter(perm_p < 0.05) %>%
  dplyr::filter(Hit > 3) %>%
  dplyr::mutate(p.value.log10 = -log10(perm_p)) %>%
  dplyr::mutate(Direction = factor(Direction, levels = c("Hyper", "Hypo")))
order <- unique(plot$Reference)
plot <- plot %>%
  dplyr::mutate(Reference = factor(Reference, levels = rev(order)))

png(file="/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/TWASEWAS_enriched_results.png", width=1600, height=1800)
ggplot(plot, aes_string(x="Actual_pct", y="Reference", size="Actual_pct", color="p.value.log10")) + 
  geom_point() +
  scale_color_continuous(low="blue", high="red", name = "-Log10(P-value)", guide=guide_colorbar(reverse=F)) +
  ylab(NULL) + ggtitle("") + DOSE::theme_dose(24) + scale_size(range=c(6, 16)) +
  facet_grid(. ~ Direction)
dev.off()





