#-----------------------------------------------------------------LIBRARIES------------------------------------------------------

library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

setwd('YOUR PROJECT PATH')

source('./scripts/helpers.R')
#-----------------------------------------------------------------COLORS------------------------------------------------------

theme_set(theme_cowplot())

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'

#-----------------------------------------------------------------READ DATA------------------------------------------------------

path_to_figures <- './figures_am/'

adrenal_medulla <- readRDS('./rds/adrenal_medulla.rds')

#-----------------------------------------------------------------SELECT DATASETS------------------------------------------------------

# Calculate number of cells coming from each batch (dataset)
cellcounts_am <- adrenal_medulla@meta.data %>% 
  group_by(orig.ident) %>% 
  summarise(Count = n())

ggplot(cellcounts_am, aes(x = orig.ident, y = Count)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Delete GSM4734604 because it has extraadrenal origin and 
# shouldn't contain cells of interest
adrenal_medulla <- subset(adrenal_medulla, subset = orig.ident != 'GSM4734604')

#-----------------------------------------------------------------QUALITY------------------------------------------------------

qcparams <- c("nFeature_RNA", "nCount_RNA", "prc.mt")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = adrenal_medulla, 
                features = qcparams[i], 
                group.by = "orig.ident", 
                pt.size = 0, raster = TRUE))
  
  ggsave2(paste0(qcparams[i], "_Vlnplot.png"), 
          path = path_to_figures, 
          width = 65, height = 15, units = "cm")
}

#-----------------------------------------------------------------NORMALIZE DATA------------------------------------------------------

adrenal_medulla[["RNA"]] <- split(adrenal_medulla[["RNA"]], 
                                  f = adrenal_medulla$orig.ident)

adrenal_medulla <- NormalizeData(adrenal_medulla, verbose = FALSE)

adrenal_medulla <- FindVariableFeatures(adrenal_medulla, verbose = FALSE)

#-----------------------------------------------------------------ACCOUNT FOR CELL CYCLE------------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

adrenal_medulla[["RNA"]] <- JoinLayers(adrenal_medulla[["RNA"]])

adrenal_medulla <- CellCycleScoring(adrenal_medulla, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, 
                                    set.ident = TRUE)

#-----------------------------------------------------------------SCALING------------------------------------------------------

# Cell cycle genes and percentage of mitochondrial genes were regressed out

adrenal_medulla <- ScaleData(adrenal_medulla, 
                             vars.to.regress = c("S.Score", "G2M.Score", "prc.mt"), 
                             features = VariableFeatures(adrenal_medulla))

#-----------------------------------------------------------------PCA------------------------------------------------------

adrenal_medulla <- RunPCA(adrenal_medulla, 
                          features = VariableFeatures(adrenal_medulla),
                          verbose = FALSE)

ElbowPlot(adrenal_medulla, ndims = 30, reduction = "pca")

ggsave2('elbow_plot.png', 
        path = path_to_figures, 
        width = 20, height = 15, units = "cm")

adrenal_medulla <- FindNeighbors(adrenal_medulla, reduction = "pca", 
                                 dims = 1:15, 
                                 verbose = FALSE)

adrenal_medulla <- FindClusters(adrenal_medulla, 
                                resolution = 0.4, 
                                verbose = FALSE)

adrenal_medulla <- RunUMAP(adrenal_medulla, reduction = "pca", 
                           dims = 1:15, verbose = FALSE)

#-----------------------------------------------------------------PLOTS------------------------------------------------------

# UMAP plot
DimPlot(adrenal_medulla, 
        group.by='seurat_clusters', 
        cols=col_vector, raster = FALSE)
ggsave2("clusters_umap.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

# Distribution of batches across obtained cluster
plot_batch_fr(adrenal_medulla, 
              cluster_column = 'seurat_clusters', 
              batch_column = 'orig.ident', 
              palette = col_vector, normalize = FALSE)

ggsave2("batch_disribution_over_clusters.png", 
        path = path_to_figures, 
        width = 18, height = 10, units = "cm")

# Distribution of batches across obtained cluster normalised to number of cells
plot_batch_fr(adrenal_medulla, 
              cluster_column = 'seurat_clusters', 
              batch_column = 'orig.ident',
              palette = col_vector, normalize = TRUE)
ggsave2("batch_disribution_over_clusters_norm.png", 
        path = path_to_figures, 
        width = 18, height = 10, units = "cm")


#-----------------------------------------------------------------MARKER GENES------------------------------------------------------

marker_genes <- c("VIM", "CXCL12", "NR1H4", "GATA4", "HAND2", "LHX2", "COL3A1", 
                  "PRRX1", "TWIST1", "TWIST2", "TBX18", "COL1A1", "COL1A2", 
                  "COL12A1", "PAX2", "LYPD1", "LHX1", "STAR", "NR5A1", "CYP17A1",
                  "CYP11A1", "CYP21A2", "SOX10", "PLP1", "FOXD3", "STMN2", 
                  "PHOX2B", "ASCL1", "ISL1", "PRPH", "ELAVL3", "ELAVL4", "CHGA", 
                  "PNMT", "TH", "MKI67", "PECAM1", "KDR", "CAVIN2", "FLT1", 
                  "EGFL7", "PRCP", "HNF4A", "AHSG", "ITIH1", "ALDOB", "VTN", 
                  "RGS10", "FCGR1A", "CD163", "AIF1",
                  "HBB", "ALAS2", "MITF", "DCT", "PMEL", "TYR")

DotPlot(adrenal_medulla, features = marker_genes, group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave2("marker_genes_clusters.pdf", 
        path = path_to_figures, 
        width = 25, height = 15, units = "cm")

#-----------------------------------------------------------------FILTERING OUT ARTEFACT CELLS------------------------------------------------------

# Several clusters consist of cells that express markers of other cell types. 
# Filter them out and rerun the analysis.

adrenal_medulla_filtered <- subset(adrenal_medulla, 
                                   subset = seurat_clusters %in% c(0, 1, 2, 4, 5, 6, 7, 9, 10))

# Check how many cells from each dataset are present in object

table(adrenal_medulla_filtered@meta.data$orig.ident)

# GSM4446537 has only 50 cells in data. 
# This might be due to the fact that this is the earliest developmental stage. 
# I would delete this dataset also.

adrenal_medulla_filtered <- subset(adrenal_medulla_filtered, 
                                   subset = orig.ident != 'GSM4446537')

#-----------------------------------------------------------------SECOND ITERATION------------------------------------------------------

adrenal_medulla_filtered[["RNA"]] <- split(adrenal_medulla_filtered[["RNA"]], 
                                           f = adrenal_medulla_filtered$orig.ident)

adrenal_medulla_filtered <- NormalizeData(adrenal_medulla_filtered,
                                          verbose = FALSE)

adrenal_medulla_filtered <- FindVariableFeatures(adrenal_medulla_filtered,
                                                 verbose = FALSE)

adrenal_medulla_filtered[["RNA"]] <- JoinLayers(adrenal_medulla_filtered[["RNA"]])

adrenal_medulla_filtered <- CellCycleScoring(adrenal_medulla_filtered, 
                                             s.features = s.genes, 
                                             g2m.features = g2m.genes, 
                                             set.ident = TRUE)

# Cell cycle genes and percentage of mitochondrial genes were regressed out
adrenal_medulla_filtered <- ScaleData(adrenal_medulla_filtered, 
                                      vars.to.regress = c("S.Score", "G2M.Score", "prc.mt"), 
                                      features = VariableFeatures(adrenal_medulla_filtered))

adrenal_medulla_filtered <- RunPCA(adrenal_medulla_filtered, 
                                   features = VariableFeatures(adrenal_medulla_filtered),
                                   verbose = FALSE)

ElbowPlot(adrenal_medulla_filtered, ndims = 30, reduction = "pca")

ggsave2('adrenal_medulla_filtered_elbow_plot.png', 
        path = path_to_figures, 
        width = 20, height = 15, units = "cm")

adrenal_medulla_filtered[["RNA"]] <- split(adrenal_medulla_filtered[["RNA"]], 
                                           f = adrenal_medulla_filtered$orig.ident)

# Integrate batches using Seurat RPCA

adrenal_medulla_filtered <- suppressMessages(
  suppressWarnings(
    IntegrateLayers(
      object = adrenal_medulla_filtered, method = RPCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.rpca", k.weight = 80
    )
  )
)

adrenal_medulla_filtered[["RNA"]] <- JoinLayers(adrenal_medulla_filtered[["RNA"]])

adrenal_medulla_filtered <- FindNeighbors(adrenal_medulla_filtered, 
                                          reduction = "integrated.rpca", 
                                          dims = 1:15,
                                          verbose = FALSE)

adrenal_medulla_filtered <- FindClusters(adrenal_medulla_filtered, 
                                         resolution = 0.3, 
                                         verbose = FALSE)

adrenal_medulla_filtered <- RunUMAP(adrenal_medulla_filtered, 
                                    reduction = "integrated.rpca", 
                                    dims = 1:15, 
                                    verbose = FALSE)

# UMAP plot
DimPlot(adrenal_medulla_filtered, group.by='seurat_clusters', 
        cols=col_vector, raster = FALSE)

ggsave2("adrenal_medulla_filtered_clusters.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

ggsave2("adrenal_medulla_filtered_clusters.pdf", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

# Marker gene expresion per cluster
DotPlot(adrenal_medulla_filtered, 
        features = marker_genes, group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Distribution of batches across obtained cluster
plot_batch_fr(adrenal_medulla_filtered, 
              cluster_column = 'seurat_clusters', 
              batch_column = 'orig.ident', 
              palette = col_vector, normalize = FALSE)

ggsave2("adrenal_medulla_filtered_batch_disribution_over_clusters.png", 
        path = path_to_figures, 
        width = 18, height = 10, units = "cm")

# Distribution of batches across obtained cluster normalised to number of cells in each batch
plot_batch_fr(adrenal_medulla_filtered, 
              cluster_column = 'seurat_clusters', 
              batch_column = 'orig.ident', 
              palette = col_vector, normalize = TRUE)

ggsave2("adrenal_medulla_filtered_batch_disribution_over_clusters_norm.png", 
        path = path_to_figures, 
        width = 18, height = 10, units = "cm")

#-----------------------------------------------------------------SUBCLUSTERING OF CLUSTER 6------------------------------------------------------

# This cluster expressed markers of 2 cell types.
# I decided to subcluster it and try to separate cell into 2 cell types

adrenal_medulla_filtered <- FindSubCluster(adrenal_medulla_filtered, cluster = 6, 
                                           graph.name = 'RNA_snn', resolution = 0.3)

DotPlot(subset(adrenal_medulla_filtered, 
               subset = seurat_clusters == 6), 
        features = marker_genes, group.by = 'sub.cluster') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave2("marker_genes_subcluster6.pdf", 
        path = path_to_figures, 
        width = 45, height = 15, units = "cm")

DimPlot(subset(adrenal_medulla_filtered, 
               subset = seurat_clusters == 6), 
        group.by='sub.cluster', 
        cols=col_vector, raster = FALSE)

ggsave2("umap_subcluster6.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

#-----------------------------------------------------------------PLOTS 2------------------------------------------------------

# Assign cell types to clusters 

celltypes_mapping <- c(
  '0' = "Sympathoblasts",
  '1' = "Sympathoblasts",
  '2' = "Chromaffin cells",
  '3' = "SCPs",
  '4' = "SCPs",
  '5' = "Sympathoblasts MKI67+",
  '6_0' = "Chromaffin cells",
  '6_1' = "Chromaffin cells",
  '6_2' = "Sympathoblasts",
  '7' = "Sympathoblasts MKI67+",
  '8' = "Sympathoblasts"
)

adrenal_medulla_filtered@meta.data$CelltypeAnnotation <- celltypes_mapping[adrenal_medulla_filtered@meta.data$sub.cluster]

DimPlot(adrenal_medulla_filtered, 
        group.by='CelltypeAnnotation', 
        cols=col_vector, raster = FALSE)

ggsave2("adrenal_medulla_filtered_annotated.png", 
        path = path_to_figures, 
        width = 20, height = 15, units = "cm", dpi = 300)

ggsave2("adrenal_medulla_filtered_annotated.pdf", 
        path = path_to_figures, 
        width = 20, height = 15, units = "cm", dpi = 300)

marker_genes_am <- c("SOX10", "PLP1", "FOXD3", "STMN2", "PHOX2B", "ASCL1", "ISL1", 
                     "PRPH", "ELAVL3", "ELAVL4", "CHGA", "PNMT", "TH", "MKI67")

DotPlot(adrenal_medulla_filtered, 
        features = marker_genes_am, 
        group.by = 'CelltypeAnnotation') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

ggsave2("marker_genes_annotated.pdf", 
        path = path_to_figures, 
        width = 15, height = 25, units = "cm")
                                               
#-----------------------------------------------------------------SAVE RDS------------------------------------------------------

saveRDS(adrenal_medulla_filtered, './rds/adrenal_medulla_preprocessed_integrated.rds')

#-----------------------------------------------------------------FIND DIFFERENTIALLY EXPRESSED GENES IN 2 SCP CLUSTERS-----------------------

adrenal_medulla_filtered_aggregated <- AggregateExpression(adrenal_medulla_filtered, 
                                                           assay='RNA', 
                                                           group.by = c("seurat_clusters"), 
                                                           return.seurat = TRUE)

p1 <- CellScatter(adrenal_medulla_filtered_aggregated, "g3", "g4") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1)

ggsave2("aggregated_SCPs.pdf", 
        path = path_to_figures, 
        width = 15, height = 15, units = "cm")



