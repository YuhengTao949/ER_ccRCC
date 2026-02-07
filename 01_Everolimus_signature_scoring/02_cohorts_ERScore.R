#### Load packages and functions ####
library(GSVA)
library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(bulkTookit)

#### Load Data and process ####
TCGA <- qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/exprSet/vsd01A.qs")
ICGC <- qs::qread("analysis/data/ccRCC_bulk_data/RECA-EU/vsd01A.qs")
E_MTAB <- qs::qread("analysis/data/ccRCC_bulk_data/E_MTAB_1980/exprSet.qs")
GSE <- qs::qread("analysis/data/ccRCC_bulk_data/GSE167573/exprSet.qs")
ERGs <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/06_ERGs.qs")
ERGs <- list(ERGs)
names(ERGs) <- "Everolimus_Response_Score"

#### ssgsea ####
res <- list(
  TCGA_KIRC = Run_ssgsea(TCGA, ERGs),
  E_MTAB_1980 = Run_ssgsea(E_MTAB, ERGs),
  ICGC = Run_ssgsea(ICGC, ERGs),
  GSE167573 = Run_ssgsea(GSE, ERGs)
)
qs::qsave(res, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ERGs_ssgsea_list.qs")

#### Add Metadata ####
metalist <- list(
  TCGA_KIRC = qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/metadata01A.qs"),
  E_MTAB_1980 = qs::qread("analysis/data/ccRCC_bulk_data/E_MTAB_1980/metadata.qs"),
  ICGC = qs::qread("analysis/data/ccRCC_bulk_data/RECA-EU/metadata01A.qs"),
  GSE167573 = qs::qread("analysis/data/ccRCC_bulk_data/GSE167573/metadata.qs")
)

meta_ssgsea_list <- lapply(names(res), function(x){
  dat1 <- metalist[[x]]
  dat2 <- res[[x]]
  dat <- cbind(dat1, dat2)
  return(dat)
})
names(meta_ssgsea_list) <- names(res)

qs::qsave(meta_ssgsea_list, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ssgsea_meta_list.qs")

#### Add exp ####
ssgsea_exp_list <- list(
  TCGA_KIRC = cbind(res[["TCGA_KIRC"]], as.data.frame(t(TCGA))),
  E_MTAB_1980 = cbind(res[["E_MTAB_1980"]], as.data.frame(t(E_MTAB))),
  ICGC = cbind(res[["ICGC"]], as.data.frame(t(ICGC))),
  GSE167573 = cbind(res[["GSE167573"]], as.data.frame(t(GSE)))
)
qs::qsave(ssgsea_exp_list, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/03_ssgsea_exp_list.qs")

#### Survival Analysis ####
# meta_ssgsea_list <- qs::qread("data/03_Everolimus_treatment_data/0303_ERGs_ssgsea/0303_02_meta_ssgsea_list.qs")
TCGA_OS <- Find_surv_cutoff(data = meta_ssgsea_list[[1]], time = "OS.time", event = "OS", variable = "Everolimus_Response_Score")
TCGA_PFI <- Find_surv_cutoff(data = meta_ssgsea_list[[1]], time = "PFI.time", event = "PFI", variable = "Everolimus_Response_Score")
TCGA_DSS <- Find_surv_cutoff(data = meta_ssgsea_list[[1]], time = "DSS.time", event = "DSS", variable = "Everolimus_Response_Score")
qs::qsave(list(TCGA_OS, TCGA_PFI, TCGA_DSS), file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/04_surv_plot_data_TCGA.qs")

##### **Survival plots** #####
p1 <- Surv_plot(TCGA_OS, group = "group", time = "OS.time", event = "OS", cutoff = "cutoff", 
                legend.title = "ERScore", title = "TCGA KIRC (OS)")
p2 <- Surv_plot(TCGA_PFI, group = "group", time = "PFI.time", event = "PFI", cutoff = "cutoff", 
                legend.title = "ERScore", title = "TCGA KIRC (PFI)")
p3 <- Surv_plot(TCGA_DSS, group = "group", time = "DSS.time", event = "DSS", cutoff = "cutoff",
                legend.title = "ERScore", title = "TCGA KIRC (DSS)")
p <- patchwork::wrap_plots(list(p1,p2,p3), ncol = 3)
print(p)
ggsave(filename = "analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/survival_analysis_ERScore_groups.pdf", p, height = 5.38, width = 13.40)

#### Tumor progress ####
data <- meta_ssgsea_list[[1]]

p_stage <- stage_boxplot(data, 
                         group_var = "ajcc_pathologic_stage", 
                         value_var = "Everolimus_Response_Score",
                         group_levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                         y_title = "Everolimus Response Score")
p_grade <- stage_boxplot(data, 
                         group_var = "tumor_grade", 
                         value_var = "Everolimus_Response_Score",
                         group_levels = c("G1", "G2", "G3", "G4"),
                         y_title = "Everolimus Response Score")
p_t <- stage_boxplot(data, 
                     group_var = "ajcc_pathologic_t", 
                     value_var = "Everolimus_Response_Score",
                     group_levels = c("T1", "T2", "T3", "T4"),
                     y_title = "Everolimus Response Score")
p <- p_grade+p_stage+p_t
print(p)
ggsave(filename = "analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/cohorts_ERScore_tumor_progress.pdf", p, width = 8.68, height = 3.86)
