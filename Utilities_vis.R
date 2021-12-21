########################
# Plot the enrichment analysis result
########################
# p_threshold: threshold of the p value
# which_p: which permutation p value should be used (nonpar or gamma)
# min_hit: minimum number of hit for each gene set
# top_n: plot top n hits
# exclude: exclude data from some databases

plot_enrichment <- function(input_dir, figure_dir, p_threshold = 0.05, which_p = "gamma", min_hit = 5, 
                            figure_width = 1600, figure_height = 800, figure_size = 6, top_n = 10, exclude = c("Tabula Muris Senis")){
  library(readr)
  library(tidyr)
  library(ggplot2)
  output_plot <- readr::read_csv(input_dir)
  
  # Select gene sets to plot
  ## Restrict to those have raw p < threshold, permutation p < threshold, hit > min_hit
  if(which_p == "gamma"){
    plot <-                # Plot those with permutation p value < 0.05 and hit > 5
      output_plot %>%
      # dplyr::filter(P_value < p_threshold) %>%
      dplyr::filter(perm_p_gamma < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::mutate(new_p = perm_p_gamma)
  }else if(which_p == "nonpar"){
    plot <-                # Plot those with permutation p value < 0.05 and hit > 5
      output_plot %>%
      # dplyr::filter(P_value < p_threshold) %>%
      dplyr::filter(perm_p_nonpar < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::mutate(new_p = perm_p_nonpar)
  }
  
  if(!is.na(exclude)){
    plot <- plot %>%
      dplyr::filter(!(Database %in% exclude))
  }
  
  # Reformat the data, and select only top_n
  plot <- plot %>%
    dplyr::arrange(P_value) %>%
    dplyr::slice_min(order_by = P_value, n = top_n)
  
  plot <-                # Plot those with permutation p value < 0.05 and hit > 5
    plot %>%
    dplyr::mutate(p.value.log10 = -log10(new_p+min(plot$new_p[plot$new_p != 0]))) %>%
    dplyr::mutate(Actual_pct = as.numeric(Actual_pct)) %>%
    dplyr::mutate(Direction = factor(Direction, levels = c("Hyper", "Hypo")))
  order <- unique(plot$Desc)
  plot <- plot %>%
    dplyr::mutate(Desc = factor(Desc, levels = rev(order)))
  
  # Plot the figure
  png(file=figure_dir, width=figure_width, height=figure_height)
  print(ggplot(plot, aes_string(x="Actual_pct", y="Desc", size="Actual_pct", color="p.value.log10")) + 
    geom_point() +
    scale_color_continuous(low="blue", high="red", name = "-Log10(P-value)", guide=guide_colorbar(reverse=F)) +
    ylab(NULL) + ggtitle("") + DOSE::theme_dose(figure_size*4) + scale_size(range=c(figure_size, figure_size*3)) +
    facet_grid(. ~ Direction) +
    theme(strip.text = element_text(size = figure_size*4, face = "bold.italic", color = "black"),
          legend.text = element_text(size = figure_size*3),
          legend.title = element_text(size = figure_size*3)))
  dev.off()
}

topN_selector <- function(input_data, p_threshold = 0.05, which_p = "gamma", min_hit = 5, top_n = 10, 
                          exclude = c("Tabula Muris Senis")){
  if(!is.na(exclude)){
    plot <- input_data %>%
      dplyr::filter(!(Database %in% exclude))
  }
  
  # Select gene sets to plot
  ## Restrict to those have raw p < threshold, permutation p < threshold, hit > min_hit
  if(which_p == "gamma"){
    select_top <- plot %>%
      dplyr::filter(perm_p_gamma < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::arrange(P_value) %>%
      dplyr::distinct(Index) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Index)
    # plot <- plot %>%
    #   dplyr::filter(Index %in% select_top)
  }else if(which_p == "nonpar"){
    select_top <- plot %>%
      dplyr::filter(perm_p_nonpar < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::arrange(P_value) %>%
      dplyr::distinct(Index) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Index)
    # plot <- plot %>%
    #   dplyr::filter(Index %in% select_top)
  }
  return(select_top)
}


plot_enrichment_heatmap <- function(input_dir = "/Users/qiyan/Dropbox/Horvath_Lab/HorvathLabCoreMembers/Qi/ToPAGE/Enrichment_Analysis_Results/EWAS_age_Ake/Nov2021", 
                                    figure_dir = "/Users/qiyan/Dropbox/Horvath_Lab/HorvathLabCoreMembers/Qi/ToPAGE/Enrichment_Analysis_Results/EWAS_age_Ake/Nov2021/Heatmap.png",
                                    p_threshold = 0.05, which_p = "gamma", min_hit = 5, 
                                    figure_width = 2800, figure_height = 1700, top_n = 10, 
                                    exclude = c("Tabula Muris Senis"), tissue_include = c("tissue",'blood','skin','brain','cortex','liver','muscle'),
                                    tissue_name = c('All','Blood','Skin','Brain','Cortex','Liver','Muscle'),
                                    cutter = NA){
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(WGCNA)
  library(gridExtra)
  library(ComplexHeatmap)
  library(tidyverse)
  
  ###########################
  # Reformat
  ###########################
  input_folder=input_dir
  TISSUE=tissue_include
  
  plot_index <- {}
  out.all2.sig <- {}
  for (i in 1:length(TISSUE)) {
    temp_file <- paste0(input_folder, "/", TISSUE[i], "/", "Enriched_TWAS_results_", TISSUE[i], ".rds")
    temp <- readRDS(temp_file)
    temp$P_value <- as.numeric(temp$P_value)
    temp_plot_index <- topN_selector(input_data = temp, p_threshold = p_threshold, which_p = which_p, min_hit = min_hit, top_n = top_n, 
                              exclude = exclude)
    plot_index <- c(plot_index, temp_plot_index)
    # long format
    temp$Tissue_ewas <- TISSUE[i]
    out.all2.sig <- rbind(out.all2.sig, temp)
  }
  
  # Need to assign an order to records
  plot_order <- data.frame("Index" = unique(plot_index), "order" = c(1:length(unique(plot_index))))
  
  out.all.sig <- out.all2.sig %>%
    dplyr::select(Index, Reference, Organism, Tissue, Cell, Type, Note, PMID, Database, Direction, Desc, Tissue_ewas, P_value, Hit, list_size) %>%
    dplyr::mutate(overlap = paste(Hit, list_size, sep = "/")) %>%
    # dplyr::select(-c(list_size)) %>%
    dplyr::arrange(Index)
  out.all.sig <- reshape2::dcast(reshape2::melt(out.all.sig, id.vars=c("Index", "Reference", "Organism", "Tissue", "Cell", "Type", "list_size",
                                                                       "Note", "PMID", "Database", "Direction", "Desc", "Tissue_ewas")), Index+Reference+Organism+Tissue+Cell+Type+list_size+Note+PMID+Database+Direction+Desc~variable+Tissue_ewas)
  
  TISSUE.short=tissue_name
  
  vars=paste0('P_value_',TISSUE)
  vars_hit=paste0('Hit_',TISSUE)
  vars.pos=paste0(vars,'.pos')
  vars.neg=paste0(vars,'.neg')
  vars_hit.pos=paste0(vars_hit,'.pos')
  vars_hit.neg=paste0(vars_hit,'.neg')
  keep.vars    =c("Index", "Reference", "Organism", "Tissue", "Cell", "Type", "list_size", "Note", "PMID", "Database", "Desc", vars, vars_hit)
  keep.vars.pos=c("Index", "Reference", "Organism", "Tissue", "Cell", "Type", "list_size", "Note", "PMID", "Database", "Desc", vars.pos, vars_hit.pos)
  out.all.sig.pos=subset(out.all.sig,Direction=='Hyper',select=keep.vars)
  names(out.all.sig.pos)=keep.vars.pos
  #
  keep.vars=c('Index',vars,vars_hit)
  keep.vars.neg=c('Index',vars.neg,vars_hit.neg)
  out.all.sig.neg=subset(out.all.sig,Direction=='Hypo',select=c('Index',vars,vars_hit))
  names(out.all.sig.neg)=keep.vars.neg
  #
  out.all.sig=merge(by='Index',out.all.sig.pos,out.all.sig.neg)
  
  #
  #very important for the variable out.sig$Order in order to align different df
  #
  out.all.sig=out.all.sig[order(out.all.sig$Index),]
  # out.sig=data.frame(Order=1:nrow(out.all.sig),out.all.sig)

  mat.p <- out.all.sig %>%
    dplyr::select(Index,vars.pos,vars.neg) %>%
    dplyr::filter(Index %in% plot_index) %>%
    dplyr::left_join(plot_order, by = "Index") %>%
    dplyr::arrange(order)
  
  mat.hit <- out.all.sig %>%                              # This is for number of hits
    dplyr::select(Index,vars_hit.pos,vars_hit.neg) %>%
    dplyr::filter(Index %in% plot_index) %>%
    dplyr::left_join(plot_order, by = "Index") %>%
    dplyr::arrange(order)
  
  {
    # May want to cut some records based on the final plot, introduce cutter
    if(!is.na(cutter)){
      mat.p <- mat.p %>%
        dplyr::filter(!(Index %in% cutter))
      mat.hit <- mat.hit %>%
        dplyr::filter(!(Index %in% cutter))
    }
  }
  
  mat.p <- mat.p %>% dplyr::select(-order)
  mat.hit <- mat.hit %>% dplyr::select(-order)
  
  mat.p <- data.frame(apply(mat.p, 2, as.numeric))  %>%
    column_to_rownames(var = "Index")
  mat.hit <- data.frame(apply(mat.hit, 2, as.numeric))  %>%
    column_to_rownames(var = "Index")
  #alternatively, you can keep the variable Order in mat.p, map.log10P then 
  #use mat.p[,-c(1)] in the complexheatmap
  # rownames(mat.p)=out.sig$Order
  mat.log10P=-log10(mat.p)
  rownames(mat.log10P)=rownames(mat.p)
  #(3) text in the heatmap
  #
  F1<-function(x,sig.cutoff0=5E-4){
    x=ifelse(x<sig.cutoff0,format(as.numeric(x),scientic=T,digit=1),'')
    
  }
  mat.text=apply(mat.p,2,F1,sig.cutoff0=0.01)
  rownames(mat.text)=rownames(mat.p)
  
  # Combine p value w/ number of hits
  mat.text_comb <- {}
  for (k in 1:nrow(mat.text)) {
    temp <- mat.text[k,]
    temp <- ifelse(temp!="", paste0(temp, " (", mat.hit[k,], ")"), "")
    mat.text_comb <- rbind(mat.text_comb, temp)
  }
  
  #
  mylayer_fun = function(j, i, x, y, width, height, fill) {
    # since grid.text can also be vectorized
    grid.text(pindex(mat.text_comb, i, j), x, y, gp = gpar(fontsize = 20, fontface='bold'))
  }
  #
  #(4) annotation
  #
  top_annotation = HeatmapAnnotation(Class = anno_block(gp =gpar(fill=  c("hotpink","skyblue")), 
                                                        labels =  c('Increased with age','Decreased with age'),
                                                        labels_gp = gpar(col = c("black","black"),
                                                                         fontsize = 28,fontface='bold'),
                                                        height=unit(20,'mm')))
  
  par(mar=c(12,30,3,3))
  col_fun = circlize::colorRamp2(c(0, max(mat.log10P)), c("white", "#FF0000"))
  col_fun(1.4)
  col_fun(5)
  col_fun = circlize::colorRamp2(c(0, 1.4,max(mat.log10P)),
                                 c("white", "#FFC9B6FF", "#FF0000"))
  #
  colnames(mat.log10P)=c(TISSUE.short,TISSUE.short)
  #
  # match.id=match(rownames(mat.p),out.all.sig$Index)
  out.sig <- out.all.sig %>%
    dplyr::filter(Index %in% rownames(mat.log10P)) %>%
    dplyr::left_join(plot_order, by = "Index") %>%
    dplyr::arrange(order)
  
  p1.mat<-Heatmap(as.matrix(mat.log10P), 
                  name = "-log10P",
                  
                  #title
                  column_title = "EWAS-TWAS enrichment",
                  column_title_gp = gpar(fontsize = 30,fontface='bold'),
                  column_names_rot = 0,
                  column_names_centered = TRUE,
                  row_title = NULL,
                  
                  #annotation
                  left_annotation = rowAnnotation(Size = anno_barplot(as.numeric(out.sig$list_size), 
                                                                      gp=gpar(border ='black', fill="grey")), width = unit(3, "cm"),
                                                  annotation_name_gp= gpar(fontsize = 22, fontface='bold')),
                  top_annotation=top_annotation,
                  
                  # Split
                  column_split=rep(1:2, each = length(TISSUE)), # one for hypermethylation and the other one for hypo,
                  row_gap = unit(.00, "mm"),
                  column_gap = unit(.00, "mm"),
                  
                  # no cluster
                  
                  cluster_columns=F,cluster_rows=F,
                  # label
                  row_names_side = c("left"),
                  row_labels =paste0(out.sig$Index,'.',out.sig$Desc),
                  row_names_max_width = max_text_width(
                    paste0(out.sig$Index,'.',out.sig$Desc), 
                    gp = gpar(fontsize = 26)),
                  
                  #col
                  col=col_fun,
                  border=T,
                  layer_fun = mylayer_fun,
                  column_names_gp = gpar(fontsize = 28,fontface='bold'),
                  row_names_gp = gpar(fontsize = 26,fontface='bold'),
                  
                  # size of heeaetmap
                  width = unit(57, "cm"),
                  height = unit(45, "cm"),
                  
                  # Legend
                  heatmap_legend_param = list(direction = "horizontal",
                                              title_gp = gpar(fontsize = 30,fontface='bold'),
                                              labels_gp = gpar(fontsize = 30,fontface='bold'),
                                              color_bar="continuous",
                                              legend_width = unit(15, "cm"))
  )
  
  #
  png(figure_dir, width = figure_width, height = figure_height)
  draw(p1.mat, heatmap_legend_side = "bottom")
  # print(p1.mat)
  dev.off()
  
  # save the output
  readr::write_csv(out.sig, file = gsub(".png", ".csv", figure_dir))
}
