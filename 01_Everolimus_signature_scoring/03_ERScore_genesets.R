#### Load packages ####
library(GSVA)
library(tidyverse)
library(bulkTookit)

#### Load data and process ####
hall <- readGMT("analysis/resource/msigdb/h.all.v2026.1.Hs.symbols.gmt")
names(hall) <- gsub("HALLMARK_", "",names(hall))
names(hall) <- gsub("_", " ", names(hall))
ERScore_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/01_ERScore_list.qs")
exp_list <- list(TCGA_KIRC = qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/exprSet/vsd01A.qs"),
                 E_MTAB_1980 = qs::qread("analysis/data/ccRCC_bulk_data/E_MTAB_1980/exprSet.qs"),
                 ICGC = qs::qread("analysis/data/ccRCC_bulk_data/RECA-EU/vsd01A.qs"),
                 GSE167573 = qs::qread("analysis/data/ccRCC_bulk_data/GSE167573/exprSet.qs"))

#### ssgsea ####
ssgsea_gsets_list <- lapply(names(exp_list), function(x){res <- Run_ssgsea(exp_list[[x]], hall)})
names(ssgsea_gsets_list) <- names(exp_list)
ERScore_gsets_score_list <- lapply(names(exp_list), function(x){
  dat <- cbind(ERScore_list[[x]], ssgsea_gsets_list[[x]])
})
names(ERScore_gsets_score_list) <- names(exp_list)
qs::qsave(ssgsea_gsets_list, file = "analysis/data/01_Everolimus_signature_scoring/03_ERScore_genesets/01_ssgsea_gsets.qs")
qs::qsave(ERScore_gsets_score_list, file = "analysis/data/01_Everolimus_signature_scoring/03_ERScore_genesets/02_ERScore_gsets.qs")

# 与 Everolimus 作用轴最直接相关
# MTORC1 SIGNALING Everolimus直接抑制mTORC1。如果ERScore与此通路正相关，可能意味着样本产生了耐药性或代偿性激活；如果负相关，则符合药物抑制预期。
# PI3K AKT MTOR SIGNALING（上游/反馈激活常见）mTOR的上游通路。肿瘤常通过上游反馈激活来逃避mTOR抑制。
# UNFOLDED PROTEIN RESPONSE（mTOR/翻译改变可牵连 ER stress）
# PROTEIN SECRETION（分泌/翻译相关）
# REACTIVE OXYGEN SPECIES PATHWAY（应激与代谢适应）
# OXIDATIVE PHOSPHORYLATION（KIRC 代谢表型差异很关键；mTOR 也会影响线粒体功能）
# GLYCOLYSIS（mTOR‑HIF 轴、缺氧代谢重编程相关） mTOR调控细胞代谢。KIRC是一种典型的代谢性癌症（Warburg效应显著），Everolimus理应抑制糖酵解。

# KIRC 进展“核心轴”：缺氧‑血管生成（建议必选）
# ANGIOGENESIS: 肾癌是多血供肿瘤，mTOR抑制剂的一个重要作用就是抗血管生成。查看ERScore是否与血管生成复燃有关。
# HYPOXIA: 缺氧信号。mTOR调控HIF1α的翻译，这在KIRC中至关重要。
# FATTY ACID METABOLISM & ADIPOGENESIS: 肾“透明”细胞癌之所以透明，是因为富含脂质。mTOR在脂质合成中起关键作用。

# 细胞周期与增殖 (代表“癌症进展”的核心)
# E2F TARGETS: 细胞周期G1/S期转换的关键，是衡量肿瘤增殖速度的黄金指标。
# G2M CHECKPOINT: 细胞分裂相关。
# MYC TARGETS V1 & MYC TARGETS V2: MYC是mTOR下游的关键效应因子，驱动肿瘤生长。
# MITOTIC SPINDLE: 有丝分裂纺锤体，代表细胞分裂活跃度。
# DNA REPAIR
# P53 PATHWAY

# 炎症反应
# INTERFERON GAMMA RESPONSE
# INTERFERON ALPHA RESPONSE
# TNFA SIGNALING VIA NFKB
# INFLAMMATORY RESPONSE
# IL6 JAK STAT3 SIGNALING（肿瘤促炎/免疫抑制轴常见）
# （可选）IL2 STAT5 SIGNALING（T 细胞活化相关）
# （可选）COMPLEMENT、COAGULATION、ALLOGRAFT REJECTION（偏“免疫/炎症/血管”综合信号）

# 肿瘤转移
# EPITHELIAL MESENCHYMAL TRANSITION 这是肿瘤转移最关键的通路。癌细胞失去上皮细胞的粘附性，获得间质细胞的迁移和侵袭能力，从而能从原发灶脱落。
# TGF BETA SIGNALING TGF-β 是诱导 EMT 的最强效细胞因子，直接驱动转移。
# WNT BETA CATENIN SIGNALING 促进肿瘤细胞的去分化、迁移能力以及维持肿瘤干细胞特性（定植必须）。
# NOTCH SIGNALING 参与细胞间通讯，调节血管生成和肿瘤干性，促进转移定植。
# HEDGEHOG SIGNALING 与肿瘤的侵袭性和干细胞维持密切相关。

#### Correlation ####
library(ComplexHeatmap)
library(circlize)

# ============================================================
# 1. 定义行分组及基因集条目
# ============================================================
pathway_groups <- list(
  "Core Mechanistic Pathway\nof mTOR Inhibition" = c(
    "MTORC1 SIGNALING",
    "PI3K AKT MTOR SIGNALING"
  ),
  "Integrative Hallmark\nPathways of ccRCC" = c(
    "ANGIOGENESIS",
    "HYPOXIA"
  ),
  "Proliferation Driver and\nCell Division Machinery" = c(
    "E2F TARGETS",
    "G2M CHECKPOINT",
    "MYC TARGETS V1",
    "MYC TARGETS V2",
    "MITOTIC SPINDLE",
    "DNA REPAIR",
    "P53 PATHWAY"
  ),
  "Inflammatory Response" = c(
    "INTERFERON GAMMA RESPONSE",
    "INTERFERON ALPHA RESPONSE",
    "TNFA SIGNALING VIA NFKB",
    "INFLAMMATORY RESPONSE",
    "IL6 JAK STAT3 SIGNALING"
  ),
  "Tumor Metastasis" = c(
    "EPITHELIAL MESENCHYMAL TRANSITION",
    "TGF BETA SIGNALING",
    "WNT BETA CATENIN SIGNALING",
    "NOTCH SIGNALING",
    "HEDGEHOG SIGNALING"
  )
)

# 有序的基因集向量
pathway_vec <- unlist(pathway_groups, use.names = FALSE)

# 构建每个基因集 -> 所属分组 的映射
group_labels <- rep(names(pathway_groups), times = sapply(pathway_groups, length))
names(group_labels) <- pathway_vec

# ============================================================
# 2. 计算相关性和显著性
# ============================================================
cohort_names <- names(ERScore_gsets_score_list)  # 四个队列

cor_mat <- matrix(NA, nrow = length(pathway_vec), ncol = length(cohort_names),
                  dimnames = list(pathway_vec, cohort_names))
pval_mat <- cor_mat  # 同维度用于存 p 值

for (cohort in cohort_names) {
  df <- ERScore_gsets_score_list[[cohort]]
  er_score <- df[, "Everolimus_Response_Score"]
  
  for (pw in pathway_vec) {
    if (pw %in% colnames(df)) {
      test <- cor.test(er_score, df[, pw], method = "spearman", exact = FALSE)
      cor_mat[pw, cohort]  <- test$estimate
      pval_mat[pw, cohort] <- test$p.value
    }
  }
}

qs::qsave(cor_mat, file = "analysis/data/01_Everolimus_signature_scoring/03_ERScore_genesets/03_cor_mat.qs")
qs::qsave(pval_mat, file = "analysis/data/01_Everolimus_signature_scoring/03_ERScore_genesets/04_cor_pval.qs")

sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                  dimnames = dimnames(pval_mat))
sig_mat[pval_mat < 0.05]  <- "*"
sig_mat[pval_mat < 0.01]  <- "**"
sig_mat[pval_mat < 0.001] <- "***"

group_colors <- c(
  "Core Mechanistic Pathway\nof mTOR Inhibition"  = "#E64B35",
  "Integrative Hallmark\nPathways of ccRCC"        = "#4DBBD5",
  "Proliferation Driver and\nCell Division Machinery" = "#00A087",
  "Inflammatory Response"                           = "#F39B7F",
  "Tumor Metastasis"                                = "#8491B4"
)

row_anno <- rowAnnotation(
  Category = group_labels[pathway_vec],
  col = list(Category = group_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Category = list(title = "Pathway Category",
                    title_gp = gpar(fontsize = 10, fontface = "bold"),
                    labels_gp = gpar(fontsize = 8))
  )
)

row_split_factor <- factor(group_labels[pathway_vec],
                           levels = names(pathway_groups))

col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                      c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B"))

ht <- Heatmap(
  cor_mat,
  name = "Spearman\nCorrelation",
  col = col_fun,

  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  row_split = row_split_factor,
  row_gap = unit(2, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,

  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_columns = FALSE,

  left_annotation = row_anno,

  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] != "") {
      grid.text(sig_mat[i, j], x, y,
                gp = gpar(fontsize = 12, fontface = "bold", col = "black"))
    }
  },
  
  rect_gp = gpar(col = "grey80", lwd = 0.8),
  
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm"),
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  
  width = unit(8, "cm"),
  column_title = "Correlation between ERScore and Hallmark Pathways",
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)

pdf("analysis/figure/01_Everolimus_signature_scoring/03_ERScore_genesets/01_cor_heatmap.pdf", width = 10, height = 10)
draw(ht, merge_legend = TRUE, padding = unit(c(10, 10, 10, 20), "mm"))
dev.off()
