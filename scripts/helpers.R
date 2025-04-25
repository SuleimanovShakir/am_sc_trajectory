library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(slingshot)
library(tradeSeq)

#-----------------------------------------------------------------PLOT DISTRIBUTION OF BATCHES ACROSS CLUSTERS------------------------------------------------------

# Custom function to plot the distribution of batches across different clusters. This function was created to check the quality of integration.

plot_batch_fr <- function(seurat_obj, cluster_column, batch_column, palette, normalize = TRUE) {
  
  data <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(cluster_column), !!sym(batch_column)) %>%
    dplyr::group_by(!!sym(cluster_column), !!sym(batch_column)) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!sym(batch_column), values_from = count, values_fill = 0)
  
  if (normalize){
    
    data_transformed <- data %>% select(!seurat_clusters) %>% mutate(across(everything(), ~ . / mean(.)))
    data_transformed$seurat_clusters <- data$seurat_clusters
    
  } else{
    
    data_transformed <- data
    
  }
  
  batch_frac <- data_transformed %>%
    rowwise() %>%
    mutate(rowSum = sum(c_across(-!!sym(cluster_column)))) %>%
    ungroup() %>%
    mutate(across(-c(!!sym(cluster_column)), ~ . / rowSum)) %>%
    dplyr::select(!rowSum)
  
  plot_data <- batch_frac %>%
    tidyr::pivot_longer(cols = -!!sym(cluster_column), names_to = batch_column, values_to = "fraction") 
  
  ggplot(plot_data, aes(x = factor(!!sym(cluster_column)), y = fraction, fill = !!sym(batch_column))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = cluster_column, y = "Fraction", fill = batch_column) +
    scale_y_continuous() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = palette)
}

#-----------------------------------------------------------------PLOT PSEUDOTIME ON UMAP------------------------------------------------------

# Custom function to plot cells in UMAP embedding with pseudotime values

plot_pseudotime_umap <- function(curves, dimred, xlim, ylim, save = TRUE, output_path = "./figures_traj/umap_pseudo.png", pseudo_range = NULL) {
  
  umap_pseudotime <- as.data.frame(dimred)
  umap_pseudotime$pseudo <- curves@assays@data@listData$pseudotime
  
  plot_df <- umap_pseudotime
  
  p <- ggplot() +
    geom_point(
      data = if (!is.null(pseudo_range)) subset(plot_df, pseudo < pseudo_range[1] | pseudo > pseudo_range[2]) else NULL,
      aes(x = umap_1, y = umap_2),
      color = "lightgrey", size = 4) +
    geom_point(
      data = if (!is.null(pseudo_range)) subset(plot_df, pseudo >= pseudo_range[1] & pseudo <= pseudo_range[2]) else plot_df,
      aes(x = umap_1, y = umap_2, color = pseudo),
      size = 4) +
    viridis::scale_color_viridis(option = "inferno", direction = 1, limits = range(plot_df$pseudo), oob = scales::squish) +
    theme_bw() +
    xlab('UMAP1') +
    ylab('UMAP2') +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2])
  
  if (save) {
    ggsave(output_path, plot = p, width = 15, height = 14, units = "cm", dpi = 300)
  } else {
    print(p)
  }
}

#-----------------------------------------------------------------PLOT GENE EXPRESSION CHANGE ALONG TRAJECTORY------------------------------------------------------

# Custom function to changes of expression for specific genes along the trajectory that was build using Slingshot

save_gene_expression_plot <- function(model, counts, gene, filename_prefix, outdir, width = 10, height = 5, units = "cm", dpi = 300) {
  
  p <- plotSmoothers(model, counts, gene = gene) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  
  ggsave(
    filename = file.path(outdir, paste0(filename_prefix, "_", gene, "_expression.png")),
    plot = p,
    width = 10,
    height = 5,
    units = "cm",
    dpi = 300
  )
}