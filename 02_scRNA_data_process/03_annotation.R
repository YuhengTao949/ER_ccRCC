#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(scTookit)

#### Load data and processe ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/02_seurat_process/01_seu_clusters.qs", nthreads = 50)
Idents(seu) <- seu$seurat_clusters
Idents(seu) <- factor(Idents(seu), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
metadata <- seu@meta.data
p <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")
print(p)

#### First level annotation ####
markerlist <- list("Immune" = "PTPRC",
                   "Epithelial" = c("EPCAM", "KRT19", "KRT18", "KRT8", "PROM1", "ALDH1A1", "CD24"),
                   "Fibroblasts" = c("COL1A1", "COL3A1", "ACAT2", "PDGFRB", "RGS5", "FGF7", "MME", "DCN", "LUM", "GSN"),
                   "Endothelial" = c("PECAM1", "VWF", "PLVAP"))

p <- MarkDotplot(markerlist, seu)
print(p)

metadata$first_level_annotation <- "unknown"
metadata$first_level_annotation[metadata$seurat_clusters %in% c("4", "0")] <- "Epithelial_Cells"
metadata$first_level_annotation[metadata$seurat_clusters == "6"] <- "Fibroblasts"
metadata$first_level_annotation[metadata$seurat_clusters %in% c("11", "7")] <- "Endothelial_Cells"
metadata$first_level_annotation[metadata$first_level_annotation == "unknown"] <- "Immune_Cells"
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "first_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
print(p)

#### Second level annotation ####
seu_immune <- subset(seu, subset = first_level_annotation == "Immune_Cells")
Idents(seu_immune) <- seu_immune$seurat_clusters
Idents(seu_immune) <- factor(Idents(seu_immune), c("1", "2", "3", "5", "8", "9", "10", "12"))

markerlist <- list("CD4_T" = c("CD4", "IL7R"),
                   "CD8_T" = c("CD8A", "CD8B"),
                   "B" = c("MS4A1", "CD19", "CD79A", "IGHG3"),
                   "NK" = c("NKG7", "GNLY", "KLRD1", "FGFBP2", "CX3CR1"),
                   "Mono" = c("S100A8", "S100A9", "LYZ", "CD14","ITGAL","ITGAX","ITGB2"),
                   "Macro" = c("APOE", "C1QA", "C1QB", "FCGR1A","CD68", "CD163","MRC1","ITGAM"),
                   "DC" = c("CD1C", "CD1E", "FCER1A", "LILRA4", "TPM2"),
                   "Mast"= c("TPSAB1", "TPSB2", "CPA3"))

p <- MarkDotplot(markerlist, seu_immune)
print(p)

metadata$second_level_annotation <- "unknown"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("4", "0")] <- "Epithelial_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "6"] <- "Fibroblasts"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("11", "7")] <- "Endothelial_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "1"] <- "CD4_T_Cells"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("9", "3")] <- "CD8_T_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "10"] <- "B_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "2"] <- "NK"
metadata$second_level_annotation[metadata$seurat_clusters == "8"] <- "Monocytes"
metadata$second_level_annotation[metadata$seurat_clusters == "5"] <- "Macrophages"
metadata$second_level_annotation[metadata$seurat_clusters == "12"] <- "Mast_Cells"

seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "second_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
print(p)

seu <- subset(seu, subset = seurat_clusters == "9", invert = T)
seu$major_celltype <- seu$second_level_annotation
p <- DimPlot(seu, reduction = "umap", group.by = "major_celltype", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
print(p)

qs::qsave(seu, file = "analysis/data/02_scRNA_data_process/03_annotation/01_seu_anno.qs", nthreads = 50)

