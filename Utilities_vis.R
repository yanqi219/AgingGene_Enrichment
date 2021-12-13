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
      dplyr::filter(Database != exclude)
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

plot_enrichment_heatmap <- function(input_dir, figure_dir, p_threshold = 0.05, which_p = "gamma", min_hit = 5, 
                                    figure_width = 1600, figure_height = 800, figure_size = 6, top_n = 10, 
                                    exclude = c("Tabula Muris Senis")){
  library(readr)
  library(tidyr)
  library(ggplot2)
  output_plot <- readr::read_csv(input_dir)
  
  if(!is.na(exclude)){
    plot <- output_plot %>%
      dplyr::filter(Database != exclude)
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
    plot <- plot %>%
      dplyr::filter(Index %in% select_top)
  }else if(which_p == "nonpar"){
    select_top <- plot %>%
      dplyr::filter(perm_p_nonpar < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::arrange(P_value) %>%
      dplyr::distinct(Index) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Index)
    plot <- plot %>%
      dplyr::filter(Index %in% select_top)
  }

  out.all.sig <- plot %>%
    dplyr::select(Index, Reference, Organism, Tissue, Cell, Type, Note, PMID, Database, Direction, Desc, P_value, Hit, list_size) %>%
    dplyr::mutate(overlap = paste(Hit, list_size, sep = "/")) %>%
    dplyr::select(-c(Hit, list_size)) %>%
    dplyr::arrange(Index)
  out.all.sig <- reshape2::dcast(reshape2::melt(out.all.sig, id.vars=c("Index", "Reference", "Organism", "Tissue", "Cell", "Type",
                                                                       "Note", "PMID", "Database", "Direction", "Desc")), Index+Reference+Organism+Tissue+Cell+Type+Note+PMID+Database+Desc~variable+Direction)
  
  
  out.all.sig=out.all.sig[order(out.all.sig$Index),]
  out.sig=data.frame(Order=1:nrow(out.all.sig),out.all.sig)
  mat.p=subset(out.sig,select=c("P_value_Hyper", "P_value_Hypo"))
  mat.p <- data.frame(apply(mat.p, 2, as.numeric))
  #alternatively, you can keep the variable Order in mat.p, map.log10P then 
  #use mat.p[,-c(1)] in the complexheatmap
  rownames(mat.p)=out.sig$Order
  mat.log10P=-log10(mat.p)

  # anno.state=subset(out.sig,select=c(Order,Index,Category,trait))
  anno.state=out.sig
  anno.state=anno.state[order(anno.state$Order),]#make sure the order of anno.state is the same with mat.p
  #
  #(2) color for category
  #
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual_col_pals$category=rownames(qual_col_pals)
  #remove light colors
  qual_col_pals=qual_col_pals[!qual_col_pals$category %in%c('Pastel1','Pastel2'),]
  #
  #blood:#F0027F
  #
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  color.df = data.frame(tissueColor=unique(col_vector))
  
  #
  anno.cat=subset(anno.state,!duplicated(Organism))
  anno.cat$row.color=color.df$tissueColor[1:nrow(anno.cat)]
  # anno.cat$row.color[anno.cat$Category=='DNAm biomarkers']='#6A3D9A'
  anno.state=merge(by='Organism',subset(anno.cat,select=c(Organism,row.color)),anno.state)
  anno.state=anno.state[order(anno.state$Order),]
  #
  anno.state.color=as.character(anno.state$row.color)
  names(anno.state.color)=as.character(anno.state$Reference)
  anno.state.color.list=list(state = anno.state.color)
  #
  anno.category=anno.state
  anno.category.color=as.character(anno.category$row.color)
  names(anno.category.color)=as.character(anno.category$Organism)
  anno.category.color.list=list(Category = anno.category.color)
  #
  #(3) text in the heatmap
  #
  F1<-function(x,sig.cutoff0=0.05){
    x=ifelse(x<sig.cutoff0,format(as.numeric(x),scientic=T,digit=1),'')
    
  }
  mat.text=apply(mat.p,2,F1,sig.cutoff0=0.05)
  rownames(mat.text)=rownames(mat.p)
  
  #
  mylayer_fun = function(j, i, x, y, width, height, fill) {
    # since grid.text can also be vectorized
    grid.text(pindex(mat.text, i, j), x, y, gp = gpar(fontsize = 18))
  }
  #
  #(4) annotation
  #
  top_annotation = HeatmapAnnotation(Class = anno_block(gp =gpar(fill=  c("hotpink","skyblue")), 
                                                        labels =  c('Hypermethylated with age','Hypomethylated with age'),
                                                        labels_gp = gpar(col = c("black","black"),
                                                                         fontsize = 26,fontface='bold')
                                                        ,height=unit(20,'mm')))
  #
  right_annotation= rowAnnotation(Category =anno.category$Organism, 
                                  # col=anno.category.color.list,
                                  annotation_label='twas category',
                                  show_annotation_name = FALSE,
                                  annotation_legend_param = list(title_gp = gpar(fontsize = 16,fontface='bold'),
                                                                 labels_gp = gpar(fontsize = 16)),show_legend=FALSE)
  #
  #(5) bar plot
  #
  anno.state=merge(by='Index',anno.state,numsig.df,all.x=T)
  anno.state=anno.state[order(anno.state$Order),]# note the order
  left_annotation= rowAnnotation(
    N.significance = anno_barplot(anno.state$index.freq,
                                  gp=gpar(border ='black',
                                          fill=as.character(anno.state.color))))
  
  #
  column_split = rep(1:2, each = length(TISSUE))# one for hypermethylation and the other one for hypo
  #
  row_split = anno.state$Category
  #=======================================================
  #png(out.summary1.png,width=20,heigh=12,unit='in',res=300)
  par(mar=c(12,30,3,3))
  col_fun = circlize::colorRamp2(c(0, max(mat.log10P)), c("white", "#FF0000"))
  col_fun(1.4)
  col_fun(5)
  col_fun = circlize::colorRamp2(c(0, 1.4,max(mat.log10P)),
                                 c("white", "#FFC9B6FF", "#FF0000"))
  #
  colnames(mat.log10P)=c(TISSUE.short,TISSUE.short)
  #
  match.id=match(rownames(mat.p),out.sig$Order)
  our.sig=out.sig[match.id,]
  p1.mat<-Heatmap(as.matrix(mat.log10P), 
                  name = "-log10P", 
                  #title
                  column_title = "EWAS-twas enrichment",
                  column_title_gp = gpar(fontsize = 30,fontface='bold'),
                  column_names_rot = 0,
                  column_names_centered = TRUE,
                  row_title = NULL,
                  #annotation
                  left_annotation = left_annotation,
                  top_annotation=top_annotation,
                  right_annotation=right_annotation,
                  column_split=column_split,
                  row_split=row_split,
                  row_gap = unit(.00, "mm"),
                  column_gap = unit(.00, "mm"),
                  #no cluster
                  cluster_columns=F,cluster_rows=F,
                  #label
                  row_names_side = c("left"),
                  row_labels =paste0(out.sig$Index,'.',out.sig$Reference),
                  row_names_max_width = max_text_width(
                    paste0(out.sig$Index,'.',out.sig$Reference), 
                    gp = gpar(fontsize = 20)),
                  
                  #col
                  col=col_fun,
                  border=T,
                  layer_fun = mylayer_fun
                  ,
                  heatmap_legend_param = list(
                    title_gp = gpar(fontsize = 20,fontface='bold'),
                    labels_gp = gpar(fontsize = 20))
                  ,
                  column_names_gp = gpar(fontsize = 18),
                  row_names_gp = gpar(fontsize = 20)
                  
                  
  )
}



