########################
# Plot the enrichment analysis result
########################
# p_threshold: threshold of the p value
# which_p: which permutation p value should be used (nonpar or gamma)
# min_hit: minimum number of hit for each gene set

plot_enrichment <- function(input_dir, figure_dir, p_threshold = 0.05, which_p = "gamma", min_hit = 5, 
                            figure_width = 1600, figure_height = 800, figure_size = 6){
  library(readr)
  library(tidyr)
  library(ggplot2)
  output_plot <- readr::read_csv(input_dir)
  
  # Select gene sets to plot
  ## Restrict to those have raw p < threshold, permutation p < threshold, and hit > min_hit
  if(which_p == "gamma"){
    plot <-                # Plot those with permutation p value < 0.05 and hit > 5
      output_plot %>%
      dplyr::filter(P_value < p_threshold) %>%
      dplyr::filter(perm_p_gamma < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::mutate(new_p = perm_p_gamma)
  }else if(which_p == "nonpar"){
    plot <-                # Plot those with permutation p value < 0.05 and hit > 5
      output_plot %>%
      dplyr::filter(P_value < p_threshold) %>%
      dplyr::filter(perm_p_nonpar < p_threshold) %>%
      dplyr::filter(as.numeric(Hit) >= min_hit) %>%
      dplyr::mutate(new_p = perm_p_nonpar)
  }
  
  # Reformat the data
  plot <-                # Plot those with permutation p value < 0.05 and hit > 5
    plot %>%
    dplyr::mutate(p.value.log10 = -log10(new_p+min(plot$new_p[plot$new_p != 0]))) %>%
    dplyr::mutate(Actual_pct = as.numeric(Actual_pct)) %>%
    dplyr::mutate(Direction = factor(Direction, levels = c("Hyper", "Hypo")))
  order <- unique(plot$Reference)
  plot <- plot %>%
    dplyr::mutate(Reference = factor(Reference, levels = rev(order)))
  
  # Plot the figure
  png(file=figure_dir, width=figure_width, height=figure_height)
  print(ggplot(plot, aes_string(x="Actual_pct", y="Reference", size="Actual_pct", color="p.value.log10")) + 
    geom_point() +
    scale_color_continuous(low="blue", high="red", name = "-Log10(P-value)", guide=guide_colorbar(reverse=F)) +
    ylab(NULL) + ggtitle("") + DOSE::theme_dose(figure_size*4) + scale_size(range=c(figure_size, figure_size*3)) +
    facet_grid(. ~ Direction) +
    theme(strip.text = element_text(size = figure_size*4, face = "bold.italic", color = "black"),
          legend.text = element_text(size = figure_size*3),
          legend.title = element_text(size = figure_size*3)))
  dev.off()
}



