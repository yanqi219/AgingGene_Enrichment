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

###########################
# Set working dir
###########################
working_dir <- '/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/DB_TWAS/'
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
    arrange(desc(abs(Meta))) %>% 
    group_by(Group) %>% slice(1:topX_cpg)
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

################## MAIN FUNCTION ##################
working_dir <- '/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/DB_TWAS/'
annotation_file_name <- "TWASDataAnnotation.csv"

output_pos <- TWASEWAS(sig_gene_list = pos, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hyper")
output_neg <- TWASEWAS(sig_gene_list = neg, ewas_species = ewas_study_species, twas_species = NA, 
                       twas_category = NA, twas_database = NA, ewas_study = input_list, topX = topX_cpg,
                       directory = working_dir, annotation_file_name = annotation_file_name, num_permutation = 100) %>%
  dplyr::mutate(Direction = "Hypo")

output_all <- rbind(output_pos, output_neg) %>%
  arrange(perm_p_nonpar)

write.table(output_all,file = "/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/TWASEWAS_enriched_results_phylo_weightAdjusted_allEWAS_lifespan.csv",sep=',',row.names = F,quote=F)
###################################################

################## DRAW PLOT ##################
output_all <- read_csv("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/TWASEWAS_enriched_results_phylo_weightAdjusted_allEWAS_lifespan.csv")
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

png(file="/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Enrichment_Analysis_Results/EWAS_maxlifespan_Caesar/TWASEWAS_enriched_results_phylo_weightAdjusted_allEWAS_lifespan.png", width=1600, height=800)
ggplot(plot, aes_string(x="Actual_pct", y="Reference", size="Actual_pct", color="p.value.log10")) + 
  geom_point() +
  scale_color_continuous(low="blue", high="red", name = "-Log10(P-value)", guide=guide_colorbar(reverse=F)) +
  ylab(NULL) + ggtitle("") + DOSE::theme_dose(24) + scale_size(range=c(6, 16)) +
  facet_grid(. ~ Direction)
dev.off()
################################################
