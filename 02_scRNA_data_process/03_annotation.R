#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(scCustomize)
library(scRNAtoolVis)
library(COSG)

#### Load data and processe ####
seu <- qs::qread("analysis/data/01_scRNA_data_process/02_seurat_process/01_seu_clusters.qs", nthreads = 50)
Idents(seu) <- seu$seurat_clusters
metadata <- seu@meta.data
p <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")
print(p)

#### First level annotation ####
markersdata <- data.frame(
  cluster = c("Immune Cells",
              rep("Epithelial Cells",7),
              rep("Fibroblasts", 10),
              rep("Endothelial Cells", 3)),
  gene = c("PTPRC", # Immune
           "EPCAM", "KRT19", "KRT18", "KRT8", "PROM1", "ALDH1A1", "CD24", # Epithelial/Cancer
           "COL1A1", "COL3A1", "ACAT2", "PDGFRB", "RGS5", "FGF7", "MME", "DCN", "LUM", "GSN", # Fibroblast
           "PECAM1", "VWF", "PLVAP") # Endothelial
)

jjDotPlot(object = seu,
          markerGene = markersdata,
          id = "seurat_clusters",
          anno = T,
          textSize = 10,
          base_size = 10,
          plot.margin = c(4,0.1,0.1,0.1))

metadata$first_level_annotation <- "unknown"
metadata$first_level_annotation[metadata$seurat_clusters %in% c("8", "0")] <- "Epithelial_Cells"
metadata$first_level_annotation[metadata$seurat_clusters == "6"] <- "Fibroblasts"
metadata$first_level_annotation[metadata$seurat_clusters %in% c("4", "10")] <- "Endothelial_Cells"
metadata$first_level_annotation[metadata$first_level_annotation == "unknown"] <- "Immune_Cells"
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "first_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
print(p)

#### Second level annotation ####
seu_immune <- subset(seu, subset = first_level_annotation == "Immune_Cells")
Idents(seu_immune) <- seu_immune$seurat_clusters

markersdata <- data.frame(
  cluster = c(rep("CD4+ T Cells", 2),
              rep("CD8+ T Cells", 2),
              rep("B Cells", 4),
              rep("NK Cells", 5),
              rep("Monocytes", 7),
              rep("Macrophages", 8),
              rep("Mast Cells", 3),
              rep("Dendritic Cells", 5)
  ),
  gene = c("CD4", "IL7R", # CD4+T
           "CD8A", "CD8B", # CD8+T
           "MS4A1", "CD19", "CD79A", "IGHG3", # B cells
           "NKG7", "GNLY", "KLRD1", "FGFBP2", "CX3CR1", # NK
           "S100A8", "S100A9", "LYZ", "CD14","ITGAL","ITGAX","ITGB2", # Mono
           "APOE", "C1QA", "C1QB", "FCGR1A","CD68", "CD163","MRC1","ITGAM", # Macro
           "TPSAB1", "TPSB2", "CPA3", # Mast
           "CD1C", "CD1E", "FCER1A", "LILRA4", "TPM2" # DC
  )
)


jjDotPlot(object = seu_immune,
          markerGene = markersdata,
          id = "seurat_clusters",
          anno = T,
          textSize = 10,
          base_size = 10,
          plot.margin = c(4,0.5,0.5,0.5))

metadata$second_level_annotation <- "unknown"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("8", "0")] <- "Epithelial_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "6"] <- "Fibroblasts"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("4", "10")] <- "Endothelial_Cells"
metadata$second_level_annotation[metadata$seurat_clusters %in% c("12", "14")] <- "B_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "13"] <- "Mast_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "5"] <- "CD8_T_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "1"] <- "CD4_T_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "3"] <- "NK_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "9"] <- "T_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "11"] <- "Dendritic_Cells"
metadata$second_level_annotation[metadata$seurat_clusters == "2"] <- "Macrophages"
metadata$second_level_annotation[metadata$seurat_clusters == "7"] <- "Monocytes"

seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "second_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
print(p)

qs::qsave(seu, file = "analysis/data/01_scRNA_data_process/03_annotation/01_seu_markers_anno.qs", nthreads = 50)

#### FindAllMarkers ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/01_scRNA_data_process/03_annotation/01_seu_markers_anno.qs", nthreads = 50)
Idents(seu) <- seu$seurat_clusters
all_markers <- FindAllMarkers(seu,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0)
all_markers$pc_d <- all_markers$pct.1 - all_markers$pct.2
all_markers$pc_r <- all_markers$pct.1/all_markers$pct.2
qs::qsave(all_markers, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/01_scRNA_data_process/03_annotation/02_FindAllMarkers_res.qs")
write.csv(all_markers, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/01_scRNA_data_process/03_annotation/02_FindAllMarkers_res.csv")

#### Enrichment ####
