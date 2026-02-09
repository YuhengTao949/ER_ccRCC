#### Load packages and functions ####
library(GSVA)
library(tidyverse)
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
qs::qsave(res, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/01_ERScore_list.qs")

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

qs::qsave(meta_ssgsea_list, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")

#### Add exp ####
ssgsea_exp_list <- list(
  TCGA_KIRC = cbind(res[["TCGA_KIRC"]], as.data.frame(t(TCGA))),
  E_MTAB_1980 = cbind(res[["E_MTAB_1980"]], as.data.frame(t(E_MTAB))),
  ICGC = cbind(res[["ICGC"]], as.data.frame(t(ICGC))),
  GSE167573 = cbind(res[["GSE167573"]], as.data.frame(t(GSE)))
)
qs::qsave(ssgsea_exp_list, file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/03_ERScore_exp_list.qs")

