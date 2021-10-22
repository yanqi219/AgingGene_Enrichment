###########################
# Technique of Prioritizing Aging-related GEnes (ToPAGE)
###########################

###########################
# The goal here is to create a pipeline for prioritizing
# aging-related genes.
###########################

###########################
# Initialize
###########################
library(futile.logger)
flog.info("Initializing ToPAGE...")
sessionInfo()
rm(list=ls())
options(StringsAsFactor=F)

###########################
# Load R Packages
###########################
flog.info("Loading R packages...")
# library(metaflow)
library(tidyverse)
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/func2_TWASEWAS.R")
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_vis.R")

###########################
# Set working dir
###########################
working_dir <- '/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/'
setwd(working_dir)
flog.info("Working directory is %s", working_dir)

###########################
# Load input CpG list
###########################
input_list = "Eutherians_phylo_weightAdjusted_allEWAS_lifespan.csv"
topX_cpg = 500                 ####### How many CpGs do we want to keep as input list

flog.info("Reading input CpG list: %s", input_list)

input <- read.csv(file = paste("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/EWASmaxlifespan_topGene_local/topCG_Caesar/", input_list, sep = ""), header = T)
input <- input %>% # Generate a column indicate whether it's hyper or hypomethylation
  dplyr::mutate(Meta = gsub("\\(.*", "", Meta)) %>%
  # dplyr::rename(Meta = "X57.PGLS.EWAS.Log.maxAgeCaesar.OrderALL.TissueALL.N215") %>%
  dplyr::mutate(Meta = as.numeric(Meta)) %>%
  dplyr::mutate(Group = ifelse(Meta > 0, "pos", "neg"))

flog.info("Number of hypermethylated CpGs: %s", sum(input$Group == "pos")) # Check the number of input probes
flog.info("Number of hypomethylated CpGs: %s", sum(input$Group == "neg"))
if(sum(input$Group == "pos") > topX_cpg | sum(input$Group == "neg") > topX_cpg){
  flog.warn("There are more than X hyper or hypomethylated probes used as input, so restrict to top X for each group")
  input <- input %>%
    dplyr::arrange(desc(abs(Meta))) %>% 
    dplyr::group_by(Group) %>% dplyr::slice(1:topX_cpg)
}

if(sum(c("CGid", "SYMBOL") %in% colnames(input)) != 2){
  flog.error("Input file should include one column named “SYMBOL” which is the gene symbol for the corresponding probe, and one column named “CGid” which is the probe ID.")
}

###########################
# Prepare input for TWASEWAS
###########################

ewas_study_species = "human" # currently it can be human, mouse, or rat; will add more later
pos <- input %>% dplyr::filter(Group == "pos")
neg <- input %>% dplyr::filter(Group == "neg")

###########################
# Start TWASEWAS
###########################

################## TWAS ##################
working_dir_2 <- paste(working_dir, "DB_TWAS/", sep = "")
annotation_file_name <- "TWASDataAnnotation.csv"

output_pos <- TWASEWAS(sig_gene_list = pos, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir_2, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hyper")
output_neg <- TWASEWAS(sig_gene_list = neg, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir_2, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hypo")
output_all <- rbind(output_pos, output_neg) %>%
  arrange(perm_p_nonpar)

write.table(output_all,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_TWAS_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)

plot_enrichment(input_dir = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_TWAS_results_phylo_weightAdjusted_allEWAS_lifespan.csv",
                figure_dir = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_TWAS_results_phylo_weightAdjusted_allEWAS_lifespan.png",
                p_threshold = 0.05, which_p = "gamma", min_hit = 5, figure_width = 2000, figure_height = 1200, figure_size = 8)
###################################################

################## Intervention ##################
working_dir_2 <- paste(working_dir, "DB_Intervention/", sep = "")
annotation_file_name <- "InterventionDataAnnotation.csv"

output_pos <- TWASEWAS(sig_gene_list = pos, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir_2, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hyper")
output_neg <- TWASEWAS(sig_gene_list = neg, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir_2, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hypo")
output_all <- rbind(output_pos, output_neg) %>%
  arrange(perm_p_nonpar)

write.table(output_all,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_intervention_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)

plot_enrichment(input_dir = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_intervention_results_phylo_weightAdjusted_allEWAS_lifespan.csv",
                figure_dir = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_intervention_results_phylo_weightAdjusted_allEWAS_lifespan.png",
                p_threshold = 0.05, which_p = "gamma", min_hit = 5, figure_width = 2000, figure_height = 1200, figure_size = 8)
###################################################

######### Use BioThing API for annotation ##########
# https://biothings.io/
# https://docs.mygene.info/en/latest/doc/query_service.html
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")

genelist <- unique(c(pos$SYMBOL, neg$SYMBOL))

geneinfo_final <- MyGeneIO_do_annotation(input_genelist = genelist, input_type = "symbol", gene_species = "human")
write.table(geneinfo_final,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Geneinfo_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)
####################################################

############### Get OMIM information ###############
# Need to require personal OMIM api key through: https://www.omim.org/api; will take only a few minutes
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")

OMIM_set_key('O9R9PTKZSCKxWM7pXPb2jg')

genelist <- unique(c(pos$SYMBOL, neg$SYMBOL))

omim_final <- {}
for (i in 1:length(genelist)) {
  gene_name_temp = genelist[i]
  temp <- OMIM_do_gene2omim(gene_name = gene_name_temp)
  omim_final <- rbind(omim_final, temp)
}
write.table(omim_final,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/OMIM_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)
####################################################

############### Get Drug-Gene interaction information ###############
# Try to get druggable genes/genes that are known drug targets based on the DGIdb. https://academic.oup.com/nar/article/49/D1/D1144/6006193?login=true
# https://www.dgidb.org/api
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")

genelist <- unique(c(pos$SYMBOL, neg$SYMBOL))

DGIdb_final <- DGIdb_gene_to_drug(gene = genelist, fda_approved_drug = T)

write.table(DGIdb_final,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/DGIdb_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)
######################################################################

############### Get PPI information ###############
# https://www.r-bloggers.com/2012/06/obtaining-a-protein-protein-interaction-network-for-a-gene-list-in-r/
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")

genelist <- unique(c(pos$SYMBOL, neg$SYMBOL))

DGIdb_final <- DGIdb_gene_to_drug(gene = genelist, fda_approved_drug = T)

write.table(DGIdb_final,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/DGIdb_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)
######################################################################

############### Combine and generate report ###############
# Combine top probes w/ gene information
Final_gene_report <- input %>%
  dplyr::arrange(-(abs(Meta))) %>%
  dplyr::mutate(meta_rank = dense_rank(-(abs(Meta)))) %>%
  dplyr::left_join(geneinfo_final, by = "SYMBOL") %>%
  dplyr::left_join(omim_final, by = "SYMBOL") %>%
  dplyr::left_join(DGIdb_final, by = "SYMBOL")

# Calculate number of enriched gene sets per probe
twas <- read_csv("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_TWAS_results_phylo_weightAdjusted_allEWAS_lifespan.csv")
intervention <- read_csv("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/Enriched_intervention_results_phylo_weightAdjusted_allEWAS_lifespan.csv")

Stat_enriched <- function(x, enrich_dat = twas){
  Num_enriched <- sum(grepl(x["SYMBOL"], enrich_dat$Hit_genes, fixed = T))
  Index_enriched <- paste(enrich_dat$Index[which(grepl(x["SYMBOL"], enrich_dat$Hit_genes, fixed = T))], collapse = ";")
  Ref_enriched <- paste(enrich_dat$Reference[which(grepl(x["SYMBOL"], enrich_dat$Hit_genes, fixed = T))], collapse = ";")
  out <- c(Num_enriched, Index_enriched, Ref_enriched)
  return(out)
}

out_twas <- data.frame(t(apply(Final_gene_report, 1, function(x) Stat_enriched(x, enrich_dat = twas))))
colnames(out_twas) <- c("num_twas_enriched", "index_twas_enriched", "ref_twas_enriched")

out_intervention <- data.frame(t(apply(Final_gene_report, 1, function(x) Stat_enriched(x, enrich_dat = intervention))))
colnames(out_intervention) <- c("num_intervention_enriched", "index_intervention_enriched", "ref_intervention_enriched")

Final_gene_report <- cbind(Final_gene_report, out_twas, out_intervention)
write_csv(Final_gene_report, file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/FinalReport_phylo_weightAdjusted_allEWAS_lifespan.csv")
###########################################################


