#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(harmony)
library(clustree)
library(ggplot2)
library(scTookit)

#### Load data and processing ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/01_data_qc/04_seu_merge.qs", nthread = 50)

#### Seurat process ####
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()
ndim <- pcs_determine(seu)
seu <- seu %>% 
  RunHarmony(reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony", dims.use = 1:ndim) %>% 
  RunUMAP(reduction = "harmony", dims = 1:ndim, reduction.name = "umap")

#### Find clusters ####
seu <- seu %>% 
  FindNeighbors(reduction = "harmony", dims = 1:ndim) %>% 
  FindClusters(resolution = seq(0.1, 1, 0.05))
p <- clustree(seu, prefix = "RNA_snn_res.")
ggsave(filename = "analysis/figure/02_scRNA_data_process/02_seurat_process/01_clustree.pdf", p, width = 15, height = 30)

seu$seurat_clusters <- seu$RNA_snn_res.0.25
qs::qsave(seu, file = "analysis/data/02_scRNA_data_process/02_seurat_process/01_seu_clusters.qs", nthread = 50)

