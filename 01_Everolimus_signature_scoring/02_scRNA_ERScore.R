#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(DESeq2)
library(GSVA)
library(ggplot2)
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype==ct]
    pseudobulk <- rowSums(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}
Run_ssgsea <- function(exp, geneset){
  exp <- as.matrix(exp)
  sp <- ssgseaParam(exprData = exp, geneSets = geneset)
  gsva_data <- gsva(sp)
  gsva_data <- gsva_data %>% t() %>% as.data.frame()
  return(gsva_data)
}

#### Load data and processing ####
seu <- qs::qread("analysis/data/01_scRNA_data_process/03_annotation/01_seu_markers_anno.qs", nthread = 50)
ERGs <- qs::qread("analysis/data/02_Everolimus_signature_scoring/01_determine_ERGs/02_ERGs.qs")
ERGs <- list(ERGs)
names(ERGs) <- "Everolimus_Response_Score"

#### Make pseudobulk matrix ####
scmeta <- seu@meta.data
counts <- seu@assays$RNA@counts

sample <- as.factor(seu$orig.ident)
names(sample) <- rownames(scmeta)
# identical(names(sample),colnames(counts))
mat.summary <- getPseudobulk(counts, celltype = sample) %>% as.data.frame()
qs::qsave(mat.summary, file = "analysis/data/02_Everolimus_signature_scoring/02_scRNA_ERScore/01_pseudobulk_counts.qs")

#### VST ####
metadata <- data.frame(sample = colnames(mat.summary))
metadata$sample <- factor(metadata$sample, levels = metadata$sample)
dds <- DESeqDataSetFromMatrix(
  countData = mat.summary, 
  colData = metadata, 
  design = ~ 1, 
  tidy = FALSE  
)

dds <- dds[rowSums(counts(dds)) > 1, ]
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "sample")
vst_matrix <- assay(vsd) %>% as.data.frame()
qs::qsave(vst_matrix, file = "analysis/data/02_Everolimus_signature_scoring/02_scRNA_ERScore/02_pseudobulk_vst.qs")

#### ssgsea ####
ERScore_mat <- Run_ssgsea(vst_matrix, ERGs)
ERScore_mat <- ERScore_mat %>% arrange(Everolimus_Response_Score) %>% rownames_to_column("sample")

median_score <- median(ERScore_mat$Everolimus_Response_Score)
ERScore_mat$group <- ifelse(ERScore_mat$Everolimus_Response_Score > median_score, "High_ERScore", "Low_ERScore")
ERScore_mat$group <- factor(ERScore_mat$group, 
                            levels = c("High_ERScore", "Low_ERScore"))
ERScore_mat$sample <- factor(ERScore_mat$sample, levels = ERScore_mat$sample)

qs::qsave(ERScore_mat, file = "analysis/data/02_Everolimus_signature_scoring/02_scRNA_ERScore/03_ERScore_mat.qs")


p <- ggplot(ERScore_mat, aes(x = sample, 
                             y = Everolimus_Response_Score,
                             color = group)) +
  geom_point(size = 2) +
  # 添加三条四分位线
  # geom_hline(yintercept = q1,             
  #            linetype = "dashed",         
  #            color = "darkgreen",         
  #            linewidth = 0.6) +           
  geom_hline(yintercept = median_score,   
             linetype = "dashed",          
             color = "black",             
             linewidth = 0.6) +           
  # geom_hline(yintercept = q3,             
  #            linetype = "dashed",         
  #            color = "darkgreen",         
  #            linewidth = 0.6) +           
  # 设置颜色映射（新增灰色中间组）
  scale_color_manual(values = c(
    "High_ERScore" = "darkred", 
    # "Mid_ERScore" = "grey60",    # 灰色数据点
    "Low_ERScore" = "darkblue"
  ), name = "ER Score Group") +
  labs(x = "Sample", 
       y = "Everolimus Resistance Score",
       title = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)
ggsave(filename = "analysis/figure/02_Everolimus_signature_scoring/01_ERScore_dotplot.pdf", p, height = 3.63, width = 4.26)

seu$ERScore_group = "unknown"
high_samples <- as.character(ERScore_mat$sample[ERScore_mat$group == "High_ERScore"])
seu$ERScore_group <- ifelse(seu$orig.ident %in% high_samples, "High_ERScore", "Low_ERScore")

p <- DimPlot(seu, pt.size = 0.05, reduction = "umap", group.by = "second_level_annotation", split.by = "ERScore_group", 
             order = F, raster = F, label = F)+ ggsci::scale_color_d3("category20")
print(p)

qs::qsave(seu, file = 'analysis/data/02_Everolimus_signature_scoring/02_scRNA_ERScore/04_seu_group.qs', nthread = 50)



