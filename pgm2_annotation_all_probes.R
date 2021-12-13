###########################
# The goal here is to create an annotation file for all the probes on mammalian array,
# so that we don't need to run the APIs everytime
###########################
library(futile.logger)
start_time <- Sys.time()

save_file_loc = "/Users/qiyan/Dropbox/Horvath_Lab/Annotation_files/Epigenomics/"
save_file_name = "Mammalian"

load("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Ortholog_Genes/MammalianMethylationConsortium/Annotations_Amin/bg_AllSpecies_Amin.RData")
annotation <- bg_AllSpecies_Amin[["human"]]; rm(bg_AllSpecies_Amin)

# ######### Use BioThing API for annotation ##########
# flog.info("Starting BioThings...")
# # https://biothings.io/
# # https://docs.mygene.info/en/latest/doc/query_service.html
# source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")
# 
# genelist <- unique(annotation$SYMBOL)
# 
# geneinfo_final <- MyGeneIO_do_annotation(input_genelist = genelist, input_type = "symbol", gene_species = "human")
# write.table(geneinfo_final,file = paste(save_file_loc, save_file_name, "/Geneinfo_results_", save_file_name, ".csv", sep = ""),sep=',',row.names = F,quote=F)
# saveRDS(geneinfo_final, file = paste(save_file_loc, save_file_name, "/Geneinfo_results_", save_file_name, ".rds", sep = ""))
# flog.info("Done BioThings")
# ####################################################
# 
# ############### Get OMIM information ###############
# flog.info("Starting OMIM...")
# # Need to require personal OMIM api key through: https://www.omim.org/api; will take only a few minutes
# source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")
# 
# OMIM_set_key('O9R9PTKZSCKxWM7pXPb2jg')
# 
# genelist <- unique(annotation$SYMBOL)
# 
# omim_final <- {}
# for (i in 1:length(genelist)) {
#   gene_name_temp = genelist[i]
#   temp <- OMIM_do_gene2omim(gene_name = gene_name_temp)
#   omim_final <- rbind(omim_final, temp)
# }
# write.table(omim_final,file = paste(save_file_loc, save_file_name, "/OMIM_results_", save_file_name, ".csv", sep = ""),sep=',',row.names = F,quote=F)
# saveRDS(omim_final, file = paste(save_file_loc, save_file_name, "/OMIM_results_", save_file_name, ".rds", sep = ""))
# flog.info("Done OMIM")
# ####################################################

############### Get Drug-Gene interaction information ###############
flog.info("Starting DGIdb...")
# Try to get druggable genes/genes that are known drug targets based on the DGIdb. https://academic.oup.com/nar/article/49/D1/D1144/6006193?login=true
# https://www.dgidb.org/api
source("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/AgingGene_Enrichment/Utilities_geneInfo.R")

genelist <- unique(annotation$SYMBOL)

DGIdb_final <- DGIdb_gene_to_drug(gene = genelist, fda_approved_drug = T)
# Check whether these drugs are aging-related based on DrugAge
drugage <- read_csv("/Users/qiyan/Dropbox/Horvath_Lab/Onging_Project/Aging_Gene_local/Reference_and_SummaryStatistics/Raw_Tables/DrugAge/DrugAge Browse.csv")
DGIdb_final$DGIdb_AgeDrug <- apply(DGIdb_final, 1, function(x) paste(toupper(stringr::str_split(x["DGIdb_drugName"], ";")[[1]]) %in% toupper(drugage$`Compound/Formulation`), collapse = ";"))

write.table(DGIdb_final,file = paste(save_file_loc, save_file_name, "/DGIdb_results_", save_file_name, ".csv", sep = ""),sep=',',row.names = F,quote=F)
saveRDS(DGIdb_final, file = paste(save_file_loc, save_file_name, "/DGIdb_results_", save_file_name, ".rds", sep = ""))
flog.info("Done DGIdb")
######################################################################

end_time <- Sys.time()
end_time - start_time