#### Load packages ####
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(forestplot)
library(grid)
library(timeROC)
set.seed(1001)

#### Load data and process ####
ERScore_meta_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")
TCGA_dat <- ERScore_meta_list[["TCGA_KIRC"]]

#### Tumor progress ####
p_stage <- stage_boxplot(TCGA_dat, 
                         group_var = "ajcc_pathologic_stage", 
                         value_var = "Everolimus_Response_Score",
                         group_levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                         y_title = "Everolimus Response Score")
p_grade <- stage_boxplot(TCGA_dat, 
                         group_var = "tumor_grade", 
                         value_var = "Everolimus_Response_Score",
                         group_levels = c("G1", "G2", "G3", "G4"),
                         y_title = "Everolimus Response Score")
p_t <- stage_boxplot(TCGA_dat, 
                     group_var = "ajcc_pathologic_t", 
                     value_var = "Everolimus_Response_Score",
                     group_levels = c("T1", "T2", "T3", "T4"),
                     y_title = "Everolimus Response Score")
p <- p_grade+p_stage+p_t
print(p)
ggsave(filename = "analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/cohorts_ERScore_tumor_progress.pdf", p, width = 8.68, height = 3.86)

#### Survival Analysis ####
# meta_ssgsea_list <- qs::qread("data/03_Everolimus_treatment_data/0303_ERGs_ssgsea/0303_02_meta_ssgsea_list.qs")
TCGA_OS <- Find_surv_cutoff(data = TCGA_dat, time = "OS.time", event = "OS", variable = "Everolimus_Response_Score")
TCGA_PFI <- Find_surv_cutoff(data = TCGA_dat, time = "PFI.time", event = "PFI", variable = "Everolimus_Response_Score")
TCGA_DSS <- Find_surv_cutoff(data = TCGA_dat, time = "DSS.time", event = "DSS", variable = "Everolimus_Response_Score")
qs::qsave(list(TCGA_OS, TCGA_PFI, TCGA_DSS), file = "analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/01_surv_plot_data_TCGA.qs")

##### **Survival plots** #####
p1 <- Surv_plot(TCGA_OS, group = "group", time = "OS.time", event = "OS", cutoff = "cutoff", 
                legend.title = "ERScore", title = "TCGA KIRC (OS)")
p2 <- Surv_plot(TCGA_PFI, group = "group", time = "PFI.time", event = "PFI", cutoff = "cutoff", 
                legend.title = "ERScore", title = "TCGA KIRC (PFI)")
p3 <- Surv_plot(TCGA_DSS, group = "group", time = "DSS.time", event = "DSS", cutoff = "cutoff",
                legend.title = "ERScore", title = "TCGA KIRC (DSS)")
p <- patchwork::wrap_plots(list(p1,p2,p3), ncol = 3)
print(p)
ggsave(filename = "analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/survival_analysis_ERScore_TCGA_KIRC.pdf", p, height = 5.38, width = 13.40)

#### Cox ####
dat_cox <- TCGA_dat %>% dplyr::select(-c("DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")) %>% filter(!is.na(OS) & !is.na(OS.time))
colnames(dat_cox) <- c("OS", "OS_time", "Gender", "Age", "Stage", "T_stage", "N_stage", "M_stage", "Tumor_grade", "ERScore")
dat_cox[is.na(dat_cox)] <- "X"
dat_cox$Stage <- gsub(" ", "_", dat_cox$Stage)
dat_cox <- dat_cox %>% 
  dplyr::mutate(Gender = factor(Gender, levels = c("male", "female")) %>% as.numeric(),
                Age = Age %>% as.numeric(),
                Stage = factor(Stage, levels = c("Stage_I", "Stage_II", "Stage_III", "Stage_IV", "X")) %>% as.numeric(),
                T_stage = factor(T_stage, levels = c("T1", "T2", "T3", "T4", "X")) %>% as.numeric(),
                N_stage = factor(N_stage, levels = c("N0", "N1", "N2", "X")) %>% as.numeric(),
                M_stage = factor(M_stage, levels = c("M0", "M1", "X")) %>% as.numeric(),
                Tumor_grade = factor(Tumor_grade, levels = c("G1", "G2", "G3", "G4", "X")) %>% as.numeric(),
                ERScore = ERScore %>% as.numeric())


##### Uni Cox ####
uni_cox <- Unicox(cox_data = dat_cox,
                  surv_event_col = "OS",
                  surv_time_col = "OS_time",
                  vars = c("Gender", "Age", "Stage", "T_stage", "N_stage", "M_stage", "Tumor_grade", "ERScore"))
pdf("analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/unicox_OS.pdf", width = 7.95, height = 4.43)
CoxForestPlot(uni_cox)
dev.off()

##### Muti Cox #####
multi_cox <- MutiCox(cox_data = dat_cox,
                     surv_event_col = "OS",
                     surv_time_col = "OS_time",
                     vars = c("Gender", "Age", "Stage", "T_stage", "N_stage", "M_stage", "Tumor_grade", "ERScore"))
pdf("analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/multicox_OS.pdf", width = 7.95, height = 4.43)
CoxForestPlot(multi_cox)
dev.off()

#### timeROC ####
pdf("analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/timeroc_os.pdf", width = 4.43, height = 4.43)
PlottimeROC(
  dat_roc    = TCGA_dat,
  time_col   = "OS.time",
  event_col  = "OS",
  var_col    = "Everolimus_Response_Score",
  main_title = "TCGA KIRC OS"
)
dev.off()

pdf("analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/timeroc_pfs.pdf", width = 4.43, height = 4.43)
PlottimeROC(
  dat_roc    = TCGA_dat,
  time_col   = "PFI.time",
  event_col  = "PFI",
  var_col    = "Everolimus_Response_Score",
  main_title = "TCGA KIRC PFI"
)
dev.off()

pdf("analysis/figure/01_Everolimus_signature_scoring/02_cohorts_ERScore/02_ERScore_TCGA_KIRC/timeroc_dfi.pdf", width = 4.43, height = 4.43)
PlottimeROC(
  dat_roc    = TCGA_dat,
  time_col   = "DSS.time",
  event_col  = "DSS",
  var_col    = "Everolimus_Response_Score",
  main_title = "TCGA KIRC DSS"
)
dev.off()
