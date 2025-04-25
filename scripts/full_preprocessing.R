#-----------------------------------------------------------------LIBRARIES------------------------------------------------------
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

setwd('YOUR PROJECT PATH')
#-----------------------------------------------------------------COLORS------------------------------------------------------

theme_set(theme_cowplot())

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'

#-----------------------------------------------------------------READ DATA------------------------------------------------------

path_to_figures <- './figures/'

filenames <- list.files('./data')

metadata <- data.frame(
    GSM = c("GSM4446535", "GSM4446536", "GSM4446537", "GSM4446538", "GSM4446539",
            "GSM4446540", "GSM4446541", "GSM4446542", "GSM4446543", "GSM4734601",
            "GSM4734602", "GSM4734603", "GSM4734604"),
    Sample = c("week8_001", "week9_063", "week6_088", "week14_123", "week12_124",
               "week8_125", "week9_005", "week11_006", "week9_007", "week8_016",
               "week9_031_paraganglia", "week12_035", "week12_036_extraadrenal"),
    stringsAsFactors = FALSE
  )

metadata$Time <- gsub('week', '', str_split_i(metadata$Sample, '_', 1))

gex_list <- list()

for (file in filenames){
  gsm <- str_split_i(file, '_', 1)
  
  data <- Read10X_h5(paste0('./data/', file))
  
  so <- CreateSeuratObject(counts = data, assay = "RNA", min.cells = 3, 
                           min.features = 200, project = gsm)
  
  gex_list[[gsm]] <- so
  
  print(paste0(gsm, ' is finished'))
}

data.so <- merge(x = gex_list[[1]], y = c(gex_list[c(2:length(gex_list))]), 
                 add.cell.ids = names(gex_list))

#saveRDS(data.so, './rds/raw_seurat.rds')

#-----------------------------------------------------------------FILTERING------------------------------------------------------

data.so <- PercentageFeatureSet(data.so, pattern = "^MT-", col.name = "prc.mt")

qcparams <- c("nFeature_RNA", "nCount_RNA", "prc.mt")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = data.so, features = qcparams[i], 
                group.by = "orig.ident", pt.size = 0, raster = TRUE))
  
  ggsave2(paste0(qcparams[i], "_Vlnplot.png"), 
          path = path_to_figures, width = 65, height = 15, units = "cm")
}

### filtering

nFeature_lower <- 200
nFeature_upper <- 7000
nCount_upper <- 70000
mt_upper <- 25

data.so_unfiltered <- data.so

data.so <- subset(data.so, subset = nCount_RNA < nCount_upper & 
                    nFeature_RNA < nFeature_upper & 
                    nFeature_RNA > nFeature_lower & prc.mt < mt_upper)

qcparams <- c("nFeature_RNA", "nCount_RNA", "prc.mt")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = data.so, features = qcparams[i], 
                group.by = "orig.ident", pt.size = 0, raster = TRUE) + 
          ylim(0, NA))
  
  ggsave2(paste0(qcparams[i], "_Vlnplot_after_filtering.png"), 
          path = path_to_figures, width = 65, height = 15, units = "cm")
}

#-----------------------------------------------------------------READ SCRUBLET RESULTS------------------------------------------------------

# Use scrublet scores to delete potential doublet cells

scrublet <- read.csv('./adata/metadata.csv', row.names = 1)

single_cells <- scrublet %>% filter(doublet_score < 0.25) %>% rownames()

data.so <- subset(data.so, cells = single_cells)

#-----------------------------------------------------------------NORMALIZE DATA------------------------------------------------------

data.so <- NormalizeData(data.so, verbose = FALSE)

data.so <- FindVariableFeatures(data.so, verbose = FALSE)

#-----------------------------------------------------------------ACCOUNT FOR CELL CYCLE------------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data.so[["RNA"]] <- JoinLayers(data.so[["RNA"]])

data.so <- CellCycleScoring(data.so, s.features = s.genes, 
                            g2m.features = g2m.genes, set.ident = TRUE)

#-----------------------------------------------------------------SCALING------------------------------------------------------

data.so <- ScaleData(data.so, 
                     vars.to.regress = c("S.Score", "G2M.Score", "prc.mt"), 
                     features = VariableFeatures(data.so))

#-----------------------------------------------------------------PCA------------------------------------------------------

data.so <- RunPCA(data.so, features = VariableFeatures(data.so), verbose = FALSE)

# Look at elbow plot to select number of principal components
ElbowPlot(data.so, ndims = 30, reduction = "pca")

ggsave2('elbow_plot.png', 
        path = path_to_figures, 
        width = 20, height = 15, units = "cm")

data.so <- FindNeighbors(data.so, 
                         reduction = "pca", 
                         dims = 1:20, 
                         verbose = FALSE)

data.so <- FindClusters(data.so, 
                        resolution = 0.4, 
                        verbose = FALSE)

data.so <- RunUMAP(data.so, 
                   reduction = "pca", 
                   dims = 1:20, 
                   verbose = FALSE)

#saveRDS(data.so, './rds/preprocessed_clustered_seurat.rds')

#-----------------------------------------------------------------PLOTS------------------------------------------------------

# UMAP by batches
DimPlot(data.so, group.by='orig.ident', cols=col_vector, raster = FALSE)

ggsave2("datasets_umap.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

# UMAP by clusters
DimPlot(data.so, group.by='seurat_clusters', cols=col_vector, raster = FALSE)

ggsave2("clusters_umap.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

#-----------------------------------------------------------------MARKER GENES------------------------------------------------------

marker_genes <- c("VIM", "CXCL12", "NR1H4", "GATA4", "HAND2", "LHX2", "COL3A1", 
                  "PRRX1", "TWIST1", "TWIST2","TBX18", "COL1A1", "COL1A2", 
                  "COL12A1", "PAX2", "LYPD1", "LHX1", "STAR", "NR5A1", "CYP17A1",
                  "CYP11A1", "CYP21A2", "SOX10", "PLP1", "FOXD3", "STMN2", 
                  "PHOX2B", "ASCL1", "ISL1", "PRPH", "ELAVL3", "ELAVL4", "CHGA", 
                  "PNMT", "TH", "PECAM1", "KDR", "CAVIN2", "FLT1", "EGFL7",
                  "PRCP", "HNF4A", "AHSG", "ITIH1", "ALDOB", "VTN", "RGS10", 
                  "FCGR1A", "CD163", "AIF1","HBB", "ALAS2", "MITF", "DCT", 
                  "PMEL", "TYR")

DotPlot(data.so, features = marker_genes, group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave2("marker_genes_clusters.pdf", 
        path = path_to_figures, 
        width = 45, height = 15, units = "cm")

#-----------------------------------------------------------------SUBCLUSTERING OF CLUSTER 9------------------------------------------------------

# Cells in this cluster expressed markers of different cell types, 
# so I decided to subcluster it and try to assign correct labels

data.so <- FindSubCluster(data.so, cluster = 9, graph.name = 'RNA_snn', 
                          resolution = 0.4)

DotPlot(subset(data.so, subset = seurat_clusters == 9), features = marker_genes, 
        group.by = 'sub.cluster') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave2("marker_genes_subcluster9.pdf", 
        path = path_to_figures, 
        width = 45, height = 15, units = "cm")

DimPlot(subset(data.so, subset = seurat_clusters == 9), 
        group.by='sub.cluster', cols=col_vector, raster = FALSE)
ggsave2("umap_subcluster9.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

FeaturePlot(data.so, features='nFeature_RNA', raster = FALSE)

#-----------------------------------------------------------------PLOTS 2------------------------------------------------------

# Assign cell types to clusters 

celltypes_mapping <- c(
  '0' = "Unknown",
  '1' = "Endothelium",
  '2' = "Adrenal gland cortex",
  '3' = "Subepicardial mesenchyme",
  '4' = "Adrenal gland cortex",
  '5' = "Erythroid cells",
  '6' = "Erythroid cells",
  '7' = "Liver primordium",
  '8' = "Sympathoblasts",
  '9_0' = "Unknown",
  '9_1' = "Unknown",
  '9_2' = "Sympathoblasts",
  '9_3' = "Chromaffin cells",
  '9_4' = "Sympathoblasts",
  '9_5' = "Unknown",
  '10' = "Endothelium",
  '11' = "Erythroid cells",
  '12' = "HSCs and immune cells",
  '13' = "Adrenal gland cortex",
  '14' = "Adrenal gland cortex",
  '15' = "Melanocytes",
  '16' = "Chromaffin cells",
  '17' = "Kidney primordium",
  '18' = "SCPs",
  '19' = "Liver primordium",
  '20' = "Intermediate mesoderm",
  '21' = "Subepicardial mesenchyme",
  '22' = "Adrenal gland cortex",
  '23' = "Endothelium",
  '24' = "HSCs and immune cells",
  '25' = "Unknown",
  '26' = "Endothelium",
  '27' = "Unknown"
)

data.so@meta.data$CelltypeAnnotation <- celltypes_mapping[data.so@meta.data$sub.cluster]

DimPlot(data.so, group.by='CelltypeAnnotation', cols=col_vector, raster = FALSE)

ggsave2("clusters_umap_annotated.png", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)

ggsave2("clusters_umap_annotated.pdf", 
        path = path_to_figures, 
        width = 30, height = 20, units = "cm", dpi = 300)


data.so@meta.data$CelltypeAnnotation <- factor(data.so@meta.data$CelltypeAnnotation, 
                                               levels = c('Intermediate mesoderm',
                                                          'Subepicardial mesenchyme',
                                                          'Kidney primordium',
                                                          'Adrenal gland cortex',
                                                          'SCPs',
                                                          'Sympathoblasts',
                                                          'Chromaffin cells',
                                                          'Endothelium',
                                                          'Liver primordium',
                                                          'HSCs and immune cells',
                                                          'Erythroid cells',
                                                          'Melanocytes',
                                                          'Unknown'))

DotPlot(data.so, features = marker_genes, group.by = 'CelltypeAnnotation') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

ggsave2("marker_genes_annotated.pdf", 
        path = path_to_figures, 
        width = 15, height = 45, units = "cm")
                                               
#-----------------------------------------------------------------SAVE RDS------------------------------------------------------

#saveRDS(data.so, './rds/preprocessed_clustered_seurat_new.rds')

# Subselect cells of interest and save to separate object
adrenal_medulla <- subset(data.so, 
                          subset = CelltypeAnnotation %in% c('SCPs', 
                                                             'Sympathoblasts', 
                                                             'Chromaffin cells'))

saveRDS(adrenal_medulla, './rds/adrenal_medulla.rds')

#-----------------------------------------------------------------FIND DIFFERENTIALLY EXPRESSED GENES PER CLUSTER--------------------------------------------

# marker_genes <- FindAllMarkers(data.so,
#                                only.pos = TRUE,
#                                min.pct = 0.4,
#                                logfc.threshold = 0.6,
#                                test.use = 'roc',
#                                group.by = 'CelltypeAnnotation')
# 
# top_markers <- marker_genes %>% 
#   group_by(cluster) %>% 
#   top_n(n = 5, wt = myAUC) %>% 
#   pull(gene)
# 
# DotPlot(data.so, features = unique(top_markers), group.by = 'CelltypeAnnotation') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   coord_flip()
# 
# ggsave2("top5_markers.pdf", 
#         path = path_to_figures, 
#         width = 15, height = 45, units = "cm")





