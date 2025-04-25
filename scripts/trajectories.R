#-----------------------------------------------------------------LIBRARIES------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(slingshot)
library(tradeSeq)

setwd('YOUR PROJECT PATH')

source('./scripts/helpers.R')

#-----------------------------------------------------------------COLORS------------------------------------------------------

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

#-----------------------------------------------------------------READ DATA------------------------------------------------------

am <- readRDS('./rds/adrenal_medulla_preprocessed_integrated.rds')

am@meta.data$CelltypeAnnotationLow <- am@meta.data$CelltypeAnnotation

# Combine 2 sympathoblasts populations for trajectories
am@meta.data$CelltypeAnnotation <- ifelse(am@meta.data$CelltypeAnnotation == 
                                            'Sympathoblasts MKI67+', 
                                          'Sympathoblasts', 
                                          am@meta.data$CelltypeAnnotation)

am_sb_ch <- subset(am, subset = CelltypeAnnotation != 'SCPs')
am_sb_scp <- subset(am, subset = CelltypeAnnotation != 'Chromaffin cells')
am_ch_scp <- subset(am, subset = CelltypeAnnotation != 'Sympathoblasts')

#-----------------------------------------------------------------EXTRACT DATA FOR TRAJECTORIES------------------------------------------------------

# extract data for:
# counts, embeddings (might be UMAP, PCA or integrated embeddings after RPCA) and cluster labels

counts_sb_ch <- as.matrix(am_sb_ch@assays$RNA$counts)
counts_sb_scp <- as.matrix(am_sb_scp@assays$RNA$counts)
counts_ch_scp <- as.matrix(am_ch_scp@assays$RNA$counts)


dimred_sb_ch <- am_sb_ch@reductions$umap@cell.embeddings
dimred_sb_scp <- am_sb_scp@reductions$umap@cell.embeddings
dimred_ch_scp <- am_ch_scp@reductions$umap@cell.embeddings


clustering_sb_ch <- am_sb_ch@meta.data$CelltypeAnnotation
clustering_sb_scp <- am_sb_scp@meta.data$CelltypeAnnotation
clustering_ch_scp <- am_ch_scp@meta.data$CelltypeAnnotation

#-----------------------------------------------------------------DEFINE LINEAGES------------------------------------------------------

# Connect cell clusters with lineage 
# These lines of code are commented because it was run on another machine

# lineages_sb_ch <- getLineages(data = dimred_sb_ch,
#                               clusterLabels = clustering_sb_ch,
#                               end.clus = "Chromaffin cells",
#                               start.clus = "Sympathoblasts")
# 
# lineages_sb_scp <- getLineages(data = dimred_sb_scp,
#                                clusterLabels = clustering_sb_scp,
#                                end.clus = "Sympathoblasts",
#                                start.clus = "SCPs")
# 
# lineages_ch_scp <- getLineages(data = dimred_ch_scp,
#                                clusterLabels = clustering_ch_scp,
#                                end.clus = "Chromaffin cells",
#                                start.clus = "SCPs")
# 
# saveRDS(lineages_sb_ch, './rds/lineages_sympathoblasts_chcells.rds')
# saveRDS(lineages_sb_scp, './rds/lineages_scps_sympathoblasts.rds')
# saveRDS(lineages_ch_scp, './rds/lineages_scps_chcells.rds')

#-----------------------------------------------------------------SMOOTH LINEAGE------------------------------------------------------

# Smooth initial lineage curve to fit the cells
# These lines of code are commented because it was run on another machine

# curves_sb_ch <- getCurves(lineages_sb_ch, approx_points = dim(lineages_sb_ch)[1])
# curves_sb_scp <- getCurves(lineages_sb_scp, approx_points = dim(lineages_sb_scp)[1])
# curves_ch_scp <- getCurves(lineages_ch_scp, approx_points = dim(lineages_ch_scp)[1])
# 
# saveRDS(curves_sb_ch, './rds/curves_sympathoblasts_chcells.rds')
# saveRDS(curves_sb_scp, './rds/curves_scps_sympathoblasts.rds')
# saveRDS(curves_ch_scp, './rds/curves_scps_chcells.rds')

#-----------------------------------------------------------------READ SAVED LINEAGES AND CURVES------------------------------------------------------

lineages_sb_ch <- readRDS('./rds/lineages_sympathoblasts_chcells.rds')
lineages_sb_scp <- readRDS('./rds/lineages_scps_sympathoblasts.rds')
lineages_ch_scp <- readRDS('./rds/lineages_scps_chcells.rds')

curves_sb_ch <- readRDS('./rds/curves_sympathoblasts_chcells.rds')
curves_sb_scp <- readRDS('./rds/curves_scps_sympathoblasts.rds')
curves_ch_scp <- readRDS('./rds/curves_scps_chcells.rds')

#-----------------------------------------------------------------PLOT PSEUDOTIME ON UMAP------------------------------------------------------

# Custom function to plot cells in UMAP embedding with pseudotime values

plot_pseudotime_umap(curves = curves_sb_ch, 
                     dimred = dimred_sb_ch, 
                     xlim = c(-10, 15), 
                     ylim = c(-15, 10), 
                     save = TRUE, 
                     output_path = './figures_traj/pseudo_sympathoblasts_chcells.png', 
                     pseudo_range = c(12, 17))

plot_pseudotime_umap(curves = curves_sb_scp, 
                     dimred = dimred_sb_scp, 
                     xlim = c(-10, 15), 
                     ylim = c(-15, 10), 
                     save = TRUE, 
                     output_path = './figures_traj/pseudo_scps_sympathoblasts.png', 
                     pseudo_range = c(4, 20))

plot_pseudotime_umap(curves = curves_ch_scp, 
                     dimred = dimred_ch_scp, 
                     xlim = c(-10, 15), 
                     ylim = c(-15, 10), 
                     save = TRUE, 
                     output_path = './figures_traj/pseudo_scps_chcells.png', 
                     pseudo_range = c(5, 19))

#-----------------------------------------------------------------GENE EXPRESSION CHANGES IN TRANSITION SYMPATHOBLASTS-CHROMAFFIN CELLS---------------

# fit generalized additive model
model_sb_ch <- fitGAM(
  counts = counts_sb_ch[VariableFeatures(am_sb_ch), ],
  sds = curves_sb_ch
)

# run the testing between pseudotime ranges
# I have fixed testing range between 12 to 17
genes_in_pseudo_range_sb_ch <- startVsEndTest(model_sb_ch, 
                                              pseudotimeValues = c(12, 17))

# filter genes by p-value and Wald statistic
genes_in_pseudo_range_sb_ch_filter <- genes_in_pseudo_range_sb_ch %>% 
  tibble::rownames_to_column("gene") %>% 
  filter(pvalue < 0.05) %>% 
  filter(abs(logFClineage1) > 1.5) %>%
  arrange(desc(waldStat))

top100_sb_ch <- head(genes_in_pseudo_range_sb_ch_filter, 100)

#-----------------------------------------------------------------GENE EXPRESSION CHANGES IN TRANSITION SCP-CHROMAFFIN CELLS---------------

# fit generalized additive model
model_scp_ch <- fitGAM(
  counts = counts_ch_scp[VariableFeatures(am_ch_scp), ],
  sds = curves_ch_scp
)

# run the testing between pseudotime ranges
# I have fixed testing range between 5 to 19
genes_in_pseudo_range_scp_ch <- startVsEndTest(model_scp_ch, 
                                               pseudotimeValues = c(5, 19))

# filter genes by p-value and Wald statistic
genes_in_pseudo_range_scp_ch_filter <- genes_in_pseudo_range_scp_ch %>% 
  tibble::rownames_to_column("gene") %>% 
  filter(pvalue < 0.05) %>% 
  filter(abs(logFClineage1) > 1.5) %>%
  arrange(desc(waldStat))

top100_scp_ch <- head(genes_in_pseudo_range_scp_ch_filter, 100)

#-----------------------------------------------------------------GENE EXPRESSION CHANGES IN TRANSITION SCP-SYMPATHOBLASTS-----------------------------

# fit generalized additive model
model_scp_sb <- fitGAM(
  counts = counts_sb_scp[VariableFeatures(am_sb_scp), ],
  sds = curves_sb_scp
)

# run the testing between pseudotime ranges
# I have fixed testing range between 4 to 20
genes_in_pseudo_range_scp_sb <- startVsEndTest(model_scp_sb, 
                                               pseudotimeValues = c(4, 20))

# filter genes by p-value and Wald statistic
genes_in_pseudo_range_scp_sb_filter <- genes_in_pseudo_range_scp_sb %>% 
  tibble::rownames_to_column("gene") %>% 
  filter(pvalue < 0.05) %>% 
  filter(abs(logFClineage1) > 1.5) %>%
  arrange(desc(waldStat))

top100_scp_sb <- head(genes_in_pseudo_range_scp_sb_filter, 100)

#-----------------------------------------------------------------UNIQUE GENES BETWEEN 2 TRANSITIONS---------------------------------------------------------

# Subselect unique genes for SCP-sympathoblasts and 
# SCP-chromaffin cells trajectories among top 100 genes from each one

genes_specific_scp_sb <- setdiff(top100_scp_sb$gene, top100_scp_ch$gene)
genes_specific_scp_ch <- setdiff(top100_scp_ch$gene, top100_scp_sb$gene)

#-----------------------------------------------------------------INTERSECTION OF GENES BETWEEN 2 TRANSITIONS---------------------------------------------------------

# Subselect genes that were identified in both transitions involving chromaffin cells

genes_specific_ch <- intersect(top100_sb_ch$gene, top100_scp_ch$gene)

#-----------------------------------------------------------------VISUALIZE SPECIFIC GENES ALONG PSEUDOTIME-----------------------------------------------

# Sympathoblasts - Chromaffin cells transition

genes_sb_ch <- c('CHGA', 'HTATSF1', 'JUNB', 'PLXNA4')

for (gene in genes_sb_ch){
  save_gene_expression_plot(model = model_sb_ch, 
                            gene = gene, 
                            counts = counts_sb_ch, 
                            filename_prefix = 'sb_ch', 
                            outdir = "./figures_traj")
  
  print(paste0('Plot for ', gene, ' is saved'))
}

# SCPs - Chromaffin cells transition

genes_scp_ch <- c('DLK1', 'PENK', 'CDKN1C', 'PLP1')

for (gene in genes_scp_ch){
  save_gene_expression_plot(model = model_scp_ch, 
                            gene = gene, 
                            counts = counts_ch_scp, 
                            filename_prefix = 'scp_ch', 
                            outdir = "./figures_traj")
  
  print(paste0('Plot for ', gene, ' is saved'))
}

# SCPs - Sympathoblasts transition

genes_scp_sb <- c('CD24', 'ELAVL4', 'SPARC', 'METRN')

for (gene in genes_scp_sb){
  save_gene_expression_plot(model = model_scp_sb, 
                            gene = gene, 
                            counts = counts_sb_scp, 
                            filename_prefix = 'scp_sb', 
                            outdir = "./figures_traj")
  
  print(paste0('Plot for ', gene, ' is saved'))
}

#-----------------------------------------------------------------VISUALIZE GENES WITH DOTPLOT---------------------------------------------------------------

genes <- c('PLP1', 'SPARC', 'METRN', 'JUNB', 'CHGA', 'CHGB', 'DLK1', 'CDKN1C', 
           'PENK', 'HTATSF1', 'MEG8', 'QDPR', 'VEGFA', 'NT5DC2', 'PLXNA4', 
           'CD24', 'ELAVL4', 'ISL1', 'STMN4', 'LINC00682', 'RTN1', 
           'MAP7', 'VSTM2L', 'MIAT', 'RBFOX1')

DotPlot(am, features = genes, group.by = 'CelltypeAnnotationLow') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab('') + 
  ylab('')

ggsave('./figures_traj/selected_genes_dotplot.pdf', 
       width = 25, height = 9, units = "cm", dpi = 300)



