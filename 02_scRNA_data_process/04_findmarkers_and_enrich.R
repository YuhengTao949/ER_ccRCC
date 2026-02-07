#### Load packages ####
library(Seurat)
library(tidyverse)

#### Load data ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/02_seurat_process/01_seu_clusters.qs", nthreads = 50)

#### FindAllMarkers ####
Idents(seu) <- seu$seurat_clusters
all_markers <- FindAllMarkers(seu,
                              only.pos = F,
                              min.pct = 0.25,
                              logfc.threshold = 0)
all_markers$pc_d <- all_markers$pct.1 - all_markers$pct.2
all_markers$pc_r <- all_markers$pct.1/all_markers$pct.2
qs::qsave(all_markers, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/04_findmarkers_and_enrich/01_findallmarkers.qs")
write.csv(all_markers, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/04_findmarkers_and_enrich/01_findallmarkers.csv")
