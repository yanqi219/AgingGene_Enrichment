####################
# Get OMIM info based on OMIM api
####################
{
  # https://www.dataquest.io/blog/r-api-tutorial/
  # https://github.com/davetang/romim
  library(XML)
  library(magrittr)
  library(xml2)
  
  OMIM_do_gene2omim <- function(gene_name = x){
    tryCatch({
      print(paste("Now searching:", gene_name))
      my_list <- OMIM_gene_to_omim(gene_name)
      my_list_omim <- sapply(my_list, OMIM_get_omim)
      temp_result <- sapply(my_list_omim, OMIM_get_title) %>% 
        data.frame(omim_id=names(.), omim_pheno=., row.names=NULL) %>%
        dplyr::mutate(omim_pheno = gsub(";", "", omim_pheno)) %>% # Remove ;
        dplyr::mutate(omim_comb = paste(omim_id, omim_pheno, sep = ":")) %>%
        dplyr::pull(omim_comb) %>%
        paste(., collapse = ";") %>%
        data.frame(SYMBOL = gene_name, OMIM_pheno = .)
      return(temp_result)},
      error=function(cond) {
        message(paste("Gene does not seem to exist in OMIM:", gene_name))
        # Choose a return value in case of error
        temp_result <- data.frame(SYMBOL = gene_name, OMIM_pheno = NA)
        return(temp_result)
      }
    )
  }
  
  OMIM_gene_to_omim <- function(gene_symbol = 'FGFR2', show_query = FALSE){
    my_search  <- paste('https://api.omim.org/api/entry/search?search=gene_symbol:', gene_symbol, sep='')
    my_include <- 'include=geneMap'
    my_query   <- paste(my_search, my_include, my_key, sep = "&")
    if (show_query){
      message(my_query)
    }
    my_result  <- read_xml(my_query) %>% xmlParse(.)
    my_list    <- xmlToList(my_result)
    my_xml     <- my_list$searchResponse$entryList$entry$geneMap$phenotypeMapList
    as.vector(unlist(sapply(my_xml, function(x) x$phenotypeMimNumber)))
  }
  
  OMIM_get_omim <- function(omim_id = 100100, #default OMIM ID
                            show_query = FALSE,
                            text = FALSE, #Includes the text field sections with the entry.
                            existflags = FALSE, #Include the 'exists' flags with the entry (clinical synopsis, allelic variant, gene map & phenotype map).
                            allelicVariantList = FALSE, #Includes the allelic variant list with the entry.
                            clinicalSynopsis = FALSE, #Include the clinical synopsis with the entry.
                            seeAlso = FALSE, #Includes the 'see also' field with the entry.
                            referenceList = FALSE, #Include the reference list with the entry.
                            geneMap = FALSE, #Include the gene map/phenotype map data with the entry.
                            externalLinks = FALSE, #Include the external links with the entry.
                            contributors = FALSE, #Includes the 'contributors' field with the entry.
                            creationDate = FALSE, #Includes the 'creation date' field with the entry.
                            editHistory = FALSE, #Includes the 'edit history' field with the entry.
                            dates = FALSE, #Include the dates data with the entry.
                            all = FALSE #Include the above data with the entry.
  ){
    #get all the arguments of the function call
    a <- as.list(match.call())
    my_mim   <- paste('mimNumber=', omim_id, sep='')
    my_link  <- 'https://api.omim.org/api/entry?'
    my_query <- paste(my_link, my_mim, my_key, sep = "&")
    #loop through all the arguments
    for (i in names(a)){
      #skip the omid_id and blank argument
      if(!i %in% '' && !i %in% 'omim_id' && !i %in% 'show_query'){
        my_include <- paste('&', 'include=', i, sep='')
        my_query <- paste(my_query, my_include, sep='')
      }
    }
    if (show_query){
      message(my_query)
    }
    read_xml(my_query) %>% xmlParse(.)
  }
  
  OMIM_get_gene <- function(my_xml){
    my_gene_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/geneSymbols")
    xmlSApply(my_gene_node, xmlValue)
  }
  
  OMIM_get_inheritance <- function(my_xml){
    my_inheritance_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeInheritance")
    xmlSApply(my_inheritance_node, xmlValue)
  }
  
  OMIM_get_pheno_key <- function(my_xml){
    my_pheno_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeMappingKey")
    xmlSApply(my_pheno_node, xmlValue)
  }
  
  OMIM_get_title <- function(my_xml){
    my_preferred_title_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/titles/preferredTitle")
    xmlSApply(my_preferred_title_node, xmlValue)
  }
  
  OMIM_set_key <- function(key){
    my_key <<- paste('apiKey=', key, sep='')
  }
}

####################
# Get Gene info based on BioThings api
####################
# https://biothings.io/
# https://docs.mygene.info/en/latest/doc/query_service.html
# https://www.dataquest.io/blog/r-api-tutorial/
{
  library(httr)
  library(jsonlite)

  MyGeneIO_gene_query <- function(gene = genelist, type = "symbol", species = "human", 
                                  field = "name,symbol,taxid,entrezgene"){
    # gene: List of genes
    # type: type of the input (symbol, entrezgene, etc.) Details can be found in https://docs.mygene.info/en/latest/doc/query_service.html
    # species: species of interest, can use common name
    # field: output information
    gene = paste(genelist,collapse = ",")
    
    params <- paste("q=", gene, "&scopes=", type, "&fields=", field, "&species=", species, sep = "")
    res <- httr::POST('http://mygene.info/v3/query', body = params, add_headers(.headers = c("Content-Type"="application/x-www-form-urlencoded")))
    data = fromJSON(rawToChar(res$content)) %>%
      dplyr::filter(!is.na(entrezgene))
    return(data)
  }
  
  my_field = "symbol,accession.genomic,accession.protein,accession.rna,AnimalQTLdb,
  biocarta,clingen.clinical_validity.disease_label,clingen.clinical_validity.online_report,
  exac.transcript,generif.text,go.BP.category,go.BP.evidence,go.BP.term,go.CC.category,
  go.CC.evidence,go.CC.term,go.MF.category,go.MF.evidence,go.MF.term,humancyc,interpro.desc,
  interpro.short_desc,kegg,map_location,mousecyc,name,netpath,other_names,pathway.biocarta.id,
  pathway.biocarta.name,pathway.humancyc.id,pathway.humancyc.name,pathway.kegg.id,pathway.kegg.name,
  pathway.mousecyc.id,pathway.mousecyc.name,pathway.netpath.id,pathway.netpath.name,pathway.pharmgkb.id,
  pathway.pharmgkb.name,pathway.pid.id,pathway.pid.name,pathway.reactome.id,pathway.reactome.name,
  pathway.smpdb.id,pathway.smpdb.name,pathway.wikipathways.id,pathway.wikipathways.name,
  pathway.yeastcyc.id,pathway.yeastcyc.name,pid,reactome,reagent.CM-LibrX-no-seq.relationship,
  reagent.CondMedia_CM_LibrAB.relationship,reagent.GNF_hs-druggable_lenti-shRNA.relationship,
  reagent.GNF_hs-druggable_plasmid-shRNA.relationship,reagent.GNF_hs-druggable_siRNA.relationship,
  reagent.GNF_hs-GPCR_IDT-siRNA.relationship,reagent.GNF_hs-oncomine_IDT-siRNA.relationship,
  reagent.GNF_hs-ORFeome1_1_reads.relationship,reagent.GNF_hs-Origene.relationship,reagent.GNF_hs-pkinase_IDT-siRNA.relationship,
  reagent.GNF_hs_LentiORF-HA-MYC.relationship,reagent.GNF_hs_LentiORF-Jred.relationship,reagent.GNF_mm+hs-MGC.relationship,
  reagent.GNF_mm+hs_RetroCDNA.relationship,reagent.GNF_mm-GIPZ_shRNA.relationship,reagent.GNF_mm-kinase_lenti-shRNA.relationship,
  reagent.GNF_mm-kinase_plasmid-shRNA.relationship,reagent.GNF_mm-TLR_lenti_shRNA.relationship,
  reagent.GNF_Qia_hs-genome_v1_siRNA.relationship,reagent.IDT_27mer_hs_ATPase_siRNAs.relationship,
  reagent.Invitrogen_IVTHSSIPKv2.relationship,reagent.MasterSecretomicsList.relationship,reagent.NIBRI_hs-Secretome_pDEST.relationship,
  reagent.NOVART_hs-genome_siRNA.relationship,reagent.Qiagen_mouse_QMIHSINHIBv1.relationship,
  reagent.Qiagen_mouse_QMIHSMIMv1.relationship,refseq.genomic,refseq.protein,refseq.rna,smpdb,summary,Vega,wikipathways,
  wikipedia.url_stub,yeastcyc"
  MyGeneIO_gene_annotation <- function(entrez = geneid, species = "human", field = my_field){
    # entrez: List of genes
    # type: type of the input (symbol, entrezgene, etc.) Details can be found in https://docs.mygene.info/en/latest/doc/query_service.html
    # species: species of interest, can use common name
    # field: output information
    gene = paste(entrez,collapse = ",")
    
    params <- paste("ids=", gene, "&fields=", field, "&species=", species, sep = "")
    res <- httr::POST('http://mygene.info/v3/gene', body = params, add_headers(.headers = c("Content-Type"="application/x-www-form-urlencoded")))
    data = fromJSON(rawToChar(res$content))
    
    # Data needs to be reformatted; a lot more informatoin can be retrieved
    ## GO terms
    bp <- data$go$BP
    bp <- data.frame(unlist(lapply(bp, function(x) paste(x$term, collapse = ";"))))
    cc <- data$go$CC
    cc <- data.frame(unlist(lapply(cc, function(x) paste(x$term, collapse = ";"))))
    mf <- data$go$MF
    mf <- data.frame(unlist(lapply(mf, function(x) paste(x$term, collapse = ";"))))
    go <- cbind(data$`_id`, bp, cc, mf); colnames(go) <- c("_id", "bp", "cc", "mf")
    ## generif
    generif <- data$generif
    generif <- data.frame(unlist(lapply(generif, function(x) paste(x$text, collapse = ";"))))
    colnames(generif) <- "generif"
    generif <- generif %>% dplyr::mutate(generif = gsub(",", "", generif))
    generif <- cbind(data$`_id`, generif); colnames(generif) <- c("_id", "generif")
    ## pathway
    kegg <- data$pathway$kegg
    kegg_id <- data.frame(unlist(lapply(kegg, function(x) paste(x$id, collapse = ";"))))
    kegg_name <- data.frame(unlist(lapply(kegg, function(x) paste(x$name, collapse = ";"))))
    kegg <- cbind(data$`_id`, kegg_id, kegg_name); colnames(kegg) <- c("_id", "kegg_id", "kegg_name")
    kegg <- kegg %>% dplyr::mutate(kegg_name = gsub(",", " ", kegg_name))
    
    reactome <- data$pathway$reactome
    reactome_id <- data.frame(unlist(lapply(reactome, function(x) paste(x$id, collapse = ";"))))
    reactome_name <- data.frame(unlist(lapply(reactome, function(x) paste(x$name, collapse = ";"))))
    reactome <- cbind(data$`_id`, reactome_id, reactome_name); colnames(reactome) <- c("_id", "reactome_id", "reactome_name")
    reactome <- reactome %>% dplyr::mutate(reactome_name = gsub(",", " ", reactome_name))
    
    wikipathways <- data$pathway$wikipathways
    wikipathways_id <- data.frame(unlist(lapply(wikipathways, function(x) paste(x$id, collapse = ";"))))
    wikipathways_name <- data.frame(unlist(lapply(wikipathways, function(x) paste(x$name, collapse = ";"))))
    wikipathways <- cbind(data$`_id`, wikipathways_id, wikipathways_name); colnames(wikipathways) <- c("_id", "wikipathways_id", "wikipathways_name")
    wikipathways <- wikipathways %>% dplyr::mutate(wikipathways_name = gsub(",", " ", wikipathways_name))
    
    pathway <- kegg %>%
      dplyr::left_join(reactome, by = "_id") %>%
      dplyr::left_join(wikipathways, by = "_id")
    
    ## Combine
    data_reformat <- data %>%
      dplyr::select(`_id`, map_location, name, symbol) %>%
      dplyr::left_join(go, by = "_id") %>%
      dplyr::left_join(generif, by = "_id") %>%
      dplyr::left_join(pathway, by = "_id")
    
    ## Make sure there's no comma
    data_reformat <- data.frame(apply(data_reformat, 2, function(x) gsub(",", " ", x)))
      
    return(data_reformat)
  }
  
  MyGeneIO_do_annotation <- function(input_genelist, input_type = "symbol", gene_species = "human"){
    gene_info <- MyGeneIO_gene_query(gene = input_genelist, type = input_type, species = gene_species, field = "name,symbol,taxid,entrezgene")
    geneid <- gene_info$entrezgene
    gene_info_final <- MyGeneIO_gene_annotation(entrez = geneid, species = gene_species) # Many info can be retrieved, for now we will focus on GO and pathway
    gene_info_final <- gene_info_final %>% dplyr::rename(SYMBOL = "symbol")
    return(gene_info_final)
  }
}

####################
# Get Drug-Gene interaction information based on DGIdb api
####################
# https://www.dgidb.org/api
# https://www.dataquest.io/blog/r-api-tutorial/
{
  library(httr)
  library(jsonlite)
  
  DGIdb_gene_to_drug <- function(gene = genelist, interaction_sources = NA, interaction_types = NA,
                                 fda_approved_drug = T, immunotherapy = FALSE, anti_neoplastic = FALSE,
                                 clinically_actionable = FALSE, druggable_genome = FALSE, drug_resistance = FALSE,
                                 gene_categories = NA, source_trust_levels = NA){
    # gene: List of genes
    # interaction_sources: A comma delimited list of source names to include in the result set. If this field is omitted, all sources will be included.
    # interaction_types: A comma delimited list of interaction types to include in the result set. If this field is omitted, all interaction types will be included.
    # fda_approved_drug: A flag denoting whether or not to limit interactions to only the ones involving fda-approved drugs. If this field is omitted, interactions for all types of drugs will be included.
    # immunotherapy: A flag denoting whether or not to limit interactions to only the ones involving immunotherapeutic drugs. If this field is omitted, interactions for all types of drugs will be included.
    # anti_neoplastic: A flag denoting whether or not to limit interactions to only the ones involving antineoplastic drugs. If this field is omitted, interactions for all types of drugs will be included.
    # clinically_actionable: A flag denoting whether or not to limit interactions to only the ones involving clinically actionable genes. If this field is omitted, interactions for all types of genes will be included.
    # druggable_genome: A flag denoting whether or not to limit interactions to only the ones involving the durggable genome. If this field is omitted, interactions for all types of genes will be included.
    # drug_resistance: A flag denoting whether or not to limit interactions to only the ones involving drug-resistant genes. If this field is omitted, interactions for all types of genes will be included.
    # gene_categories: A comma delimited list of gene categories to include in the result set. If this field is omitted, all gene categories will be included.
    # source_trust_levels: A comma delimited list of source trust levels to include in the result set. If this field is omitted, all trust levels will be included.
    
    gene = paste(genelist,collapse = ",")
    my_query <- paste("https://dgidb.org/api/v2/interactions.json?", "genes=", gene, sep = "")
    
    # Apply filters
    if(!all(is.na(interaction_sources))){
      my_query <- paste(my_query, "&interaction_sources=", interaction_sources, sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(interaction_types))){
      my_query <- paste(my_query, "&interaction_types=", interaction_types, sep = "") %>% gsub(" ", "", .)
    }
    if(fda_approved_drug){
      my_query <- paste(my_query, "&fda_approved_drug=true", sep = "") %>% gsub(" ", "", .)
    }
    if(immunotherapy){
      my_query <- paste(my_query, "&immunotherapy=true", sep = "") %>% gsub(" ", "", .)
    }
    if(anti_neoplastic){
      my_query <- paste(my_query, "&anti_neoplastic=true", sep = "") %>% gsub(" ", "", .)
    }
    if(clinically_actionable){
      my_query <- paste(my_query, "&clinically_actionable=true", sep = "") %>% gsub(" ", "", .)
    }
    if(druggable_genome){
      my_query <- paste(my_query, "&druggable_genome=true", sep = "") %>% gsub(" ", "", .)
    }
    if(drug_resistance){
      my_query <- paste(my_query, "&drug_resistance=true", sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(gene_categories))){
      my_query <- paste(my_query, "&gene_categories=", gene_categories, sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(source_trust_levels))){
      my_query <- paste(my_query, "&source_trust_levels=", source_trust_levels, sep = "") %>% gsub(" ", "", .)
    }
    
    # Use API
    res <- httr::GET(my_query)
    data = fromJSON(rawToChar(res$content))
    
    # Obtain gene category
    geneCategories <- data$matchedTerms$geneCategories
    DGIdb_geneCategories <- data.frame(unlist(lapply(geneCategories, function(x) paste(x$name, collapse = ";"))))
    DGIdb_geneCategories_comb <- cbind(data$matchedTerms$geneName, DGIdb_geneCategories)
    colnames(DGIdb_geneCategories_comb) <- c("SYMBOL", "DGIdb_geneCategories")
    # Obtain drug interaction
    interactions <- data$matchedTerms$interactions # from interaction: keep interaction type, drug name, drug ID, sources, pmids, and score
    DGIdb_interactionTypes <- data.frame(unlist(lapply(interactions, function(x) paste(x$interactionTypes, collapse = ";"))))
    DGIdb_drugName <- data.frame(unlist(lapply(interactions, function(x) paste(x$drugName, collapse = ";"))))
    DGIdb_drugConceptId <- data.frame(unlist(lapply(interactions, function(x) paste(x$drugConceptId, collapse = ";"))))
    DGIdb_sources <- data.frame(unlist(lapply(interactions, function(x) paste(x$sources, collapse = ";"))))
    DGIdb_pmids <- data.frame(unlist(lapply(interactions, function(x) paste(x$pmids, collapse = ";"))))
    DGIdb_score <- data.frame(unlist(lapply(interactions, function(x) paste(x$score, collapse = ";"))))
    DGIdb_interaction_comb <- cbind(data$matchedTerms$geneName, DGIdb_interactionTypes, DGIdb_drugName, DGIdb_drugConceptId, DGIdb_sources,
                                    DGIdb_pmids, DGIdb_score)
    colnames(DGIdb_interaction_comb) <- c("SYMBOL", "DGIdb_interactionTypes", "DGIdb_drugName", "DGIdb_drugConceptId", "DGIdb_sources", 
                                          "DGIdb_pmids", "DGIdb_score")
    # Put together
    DGIdb_final <- DGIdb_geneCategories_comb %>% dplyr::left_join(DGIdb_interaction_comb, by = "SYMBOL")
    
    return(DGIdb_final)
  }
  
  DGIdb_drug_to_gene <- function(drug = druglist, interaction_sources = NA, interaction_types = NA,
                                 fda_approved_drug = FALSE, immunotherapy = FALSE, anti_neoplastic = FALSE,
                                 clinically_actionable = FALSE, druggable_genome = FALSE, drug_resistance = FALSE,
                                 gene_categories = NA, source_trust_levels = NA){
    # drug: List of drugs
    
    drug = paste(druglist,collapse = ",")
    my_query <- paste("https://dgidb.org/api/v2/interactions.json?", "drugs=", drug, sep = "")
    
    # Apply filters
    if(!all(is.na(interaction_sources))){
      my_query <- paste(my_query, "&interaction_sources=", interaction_sources, sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(interaction_types))){
      my_query <- paste(my_query, "&interaction_types=", interaction_types, sep = "") %>% gsub(" ", "", .)
    }
    if(fda_approved_drug){
      my_query <- paste(my_query, "&fda_approved_drug=true", sep = "") %>% gsub(" ", "", .)
    }
    if(immunotherapy){
      my_query <- paste(my_query, "&immunotherapy=true", sep = "") %>% gsub(" ", "", .)
    }
    if(anti_neoplastic){
      my_query <- paste(my_query, "&anti_neoplastic=true", sep = "") %>% gsub(" ", "", .)
    }
    if(clinically_actionable){
      my_query <- paste(my_query, "&clinically_actionable=true", sep = "") %>% gsub(" ", "", .)
    }
    if(druggable_genome){
      my_query <- paste(my_query, "&druggable_genome=true", sep = "") %>% gsub(" ", "", .)
    }
    if(drug_resistance){
      my_query <- paste(my_query, "&drug_resistance=true", sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(gene_categories))){
      my_query <- paste(my_query, "&gene_categories=", gene_categories, sep = "") %>% gsub(" ", "", .)
    }
    if(!all(is.na(source_trust_levels))){
      my_query <- paste(my_query, "&source_trust_levels=", source_trust_levels, sep = "") %>% gsub(" ", "", .)
    }
    
    # Use API
    res <- httr::GET(my_query)
    data = fromJSON(rawToChar(res$content))
    
    # Obtain drug interaction
    interactions <- data$matchedTerms$interactions # from interaction: keep interaction type, gene name, entrez ID, sources, pmids, and score
    DGIdb_interactionTypes <- data.frame(unlist(lapply(interactions, function(x) paste(x$interactionTypes, collapse = ";"))))
    DGIdb_geneName <- data.frame(unlist(lapply(interactions, function(x) paste(x$geneName, collapse = ";"))))
    DGIdb_geneEntrezId <- data.frame(unlist(lapply(interactions, function(x) paste(x$geneEntrezId, collapse = ";"))))
    DGIdb_sources <- data.frame(unlist(lapply(interactions, function(x) paste(x$sources, collapse = ";"))))
    DGIdb_pmids <- data.frame(unlist(lapply(interactions, function(x) paste(x$pmids, collapse = ";"))))
    DGIdb_score <- data.frame(unlist(lapply(interactions, function(x) paste(x$score, collapse = ";"))))
    DGIdb_interaction_comb <- cbind(data$matchedTerms$drugName, data$matchedTerms$conceptId, DGIdb_interactionTypes, DGIdb_geneName, 
                                    DGIdb_geneEntrezId, DGIdb_sources, DGIdb_pmids, DGIdb_score)
    colnames(DGIdb_interaction_comb) <- c("SYMBOL", "DGIdb_drugID", "DGIdb_interactionTypes", "DGIdb_geneName", "DGIdb_geneEntrezId", 
                                          "DGIdb_sources", "DGIdb_pmids", "DGIdb_score")
    
    return(DGIdb_interaction_comb)
  }
}

####################
# Get protein-protein interaction network based on NCBI - BioGRID and HPRD
####################
# https://www.r-bloggers.com/2012/06/obtaining-a-protein-protein-interaction-network-for-a-gene-list-in-r/
{
  PPI_getppiNCBI <- function(g.n) {
    require(XML)
    ppi <- data.frame()
    for(i in 1:length(g.n)){
      o <- htmlParse(paste("https://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
      # check if interaction table exists
      exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
      if(exist){
        p <- getNodeSet(o, "//table")
        ## need to know which table is the good one
        for(j in 1:length(p)){
          int <- readHTMLTable(p[[j]])
          if(colnames(int)[2]=="Interactant"){break}
        }
        ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
      }
      # play nice! and avoid being kicked out from NCBI servers
      Sys.sleep(1)
    }
    if(dim(ppi)[1]>0){
      ppi <- unique(ppi)
      print(paste(dim(ppi)[1], "interactions found"))
      return(ppi)
    } else{
      print("No interaction found")
    }
  }
}
