#### Environment set ####
Sys.setenv(http_proxy = "http://iyun70.com:7890")
Sys.setenv(https_proxy = "http://iyun70.com:7890")
options(timeout = 300)
#### Load packages ####
library(tidyverse)
library(DESeq2)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(rrvgo)
library(ggplot2)
library(ggrepel)
library(EnrichExplainer)

#### Load data and process ####
counts <- qs::qread("RNA_seq_upstream/GSE99875/02_annotated_expression_matrices/counts.qs") 
metadata <- qs::qread("RNA_seq_upstream/GSE99875/resource/metadata.qs")
path_hierarchy <- jsonlite::fromJSON("analysis/resource/br08901.json", simplifyDataFrame = T, flatten = T)
kegg_df <- path_hierarchy$children %>% 
  rename(level1 = name) %>%
  unnest(cols = c(children)) %>% 
  rename(level2=name) %>%
  unnest(cols = c("children")) %>%
  rename(level3=name) %>%
  mutate(id = paste0("hsa", substr(level3, 1,5)),
         level3 = substr(level3, 7, nchar(level3))) %>%
  select(id, level3, level2, level1)

#### Identify differentially expressed genes ####
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group)
dds <- DESeq(dds)
contrast = c("group", "Everolimus", "Control")
dd1 <- results(dds, contrast = contrast, alpha = 0.05)
plotMA(dd1, ylim = c(-10, 10))
dd2 <- lfcShrink(dds, contrast = contrast, res = dd1, type = "ashr")
plotMA(dd2, ylim = c(-10, 10))
deg_res <- dd2 %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  rownames_to_column("gene")
qs::qsave(deg_res,file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/01_DESeq2_res.qs")

#### Prepare for LLM ####
deg_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/01_DESeq2_res.qs")
deg_filter <- deg_res %>% dplyr::filter(padj < 0.05 & log2FoldChange > 1)
diff_genes <- deg_filter$log2FoldChange
names(diff_genes) <- deg_filter$gene
diff_genes <- sort(diff_genes, decreasing = T)

#### Enrichment ####
gene_mapping <- bitr(deg_res$gene, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = "org.Hs.eg.db")
colnames(gene_mapping) <- c("gene", "entrez")
gene_df <- deg_res %>% dplyr::inner_join(gene_mapping, by = "gene")

univers_entrez <- gene_df %>% dplyr::filter(baseMean > 0) %>% pull(entrez)
entrez <- gene_df %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange > 1) %>% 
  pull(entrez)

entrez_list <- gene_df$log2FoldChange
names(entrez_list) <- gene_df$entrez
entrez_list <- sort(entrez_list, decreasing = TRUE)

##### GO #####
ora_go_res <- enrichGO(gene          = entrez,
                       keyType       = "ENTREZID",
                       universe      = univers_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       minGSSize     = 0,
                       maxGSSize     = 1000000,
                       pvalueCutoff  = 0.05)
ora_go_res <- setReadable(ora_go_res,
                          OrgDb = "org.Hs.eg.db",
                          keyType = "ENTREZID")
ora_go_res@result <- subset(ora_go_res@result, p.adjust < 0.05)
ora_go_dat <- as.data.frame(ora_go_res) %>% dplyr::arrange(desc(FoldEnrichment))
qs::qsave(ora_go_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/02_ora_go.qs")

gsea_go_res <- gseGO(geneList     = entrez_list,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     minGSSize    = 0,
                     maxGSSize    = 1000000,
                     pvalueCutoff = 0.05,
                     verbose      = TRUE)
gsea_go_res <- setReadable(gsea_go_res,
                           OrgDb = "org.Hs.eg.db",
                           keyType = "ENTREZID")
gsea_go_res@result <- subset(gsea_go_res@result, NES > 0 & p.adjust < 0.05)
gsea_go_dat <- as.data.frame(gsea_go_res) %>% arrange(desc(NES))
qs::qsave(gsea_go_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/02_gsea_go.qs")

##### rrvgo #####
simMatrix <- calculateSimMatrix(gsea_go_dat$ID,
                                orgdb = "org.Hs.eg.db",
                                ont = "BP",
                                method = "Wang")
scores <- setNames(gsea_go_dat$NES, gsea_go_dat$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
qs::qsave(reducedTerms, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/02_rrvgo.qs")

gsea_go_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/02_gsea_go.qs")
# reducedTerms <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/02_rrvgo.qs")
# gsea_go_res@result <- subset(gsea_go_res@result, ID %in% reducedTerms$parent)
prompt <- interpret_tool(gsea_go_res, database = "GO", pathway_num = 50, diff_genes, contact_LLM = FALSE,
                         context = "786-0 cells were treated with either a control DMSO vehicle or 10um Everolimus, an mTOR inhibitor, for 4 hours.")
writeLines(prompt, "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/prompt/gsea_go.txt")

##### KEGG #####
ora_kegg_res <- enrichKEGG(gene         = entrez,
                           organism     = 'hsa',
                           pvalueCutoff = 0.05,
                           universe     = univers_entrez,
                           minGSSize    = 0,
                           maxGSSize    = 1000000)
ora_kegg_res <- setReadable(ora_kegg_res,
                            OrgDb = "org.Hs.eg.db",
                            keyType = "ENTREZID")
ora_kegg_res@result <- subset(ora_kegg_res@result, p.adjust < 0.05)
ora_kegg_dat <- as.data.frame(ora_kegg_res) %>% dplyr::arrange(desc(FoldEnrichment))
qs::qsave(ora_kegg_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/03_ora_kegg.qs")

gsea_kegg_res <- gseKEGG(geneList     = entrez_list,
                         organism     = 'hsa',
                         minGSSize    = 0,
                         maxGSSize    = 1000000,
                         pvalueCutoff = 0.05,
                         verbose      = TRUE)
gsea_kegg_res <- setReadable(gsea_kegg_res,
                             OrgDb = "org.Hs.eg.db",
                             keyType = "ENTREZID")
gsea_kegg_res@result <- subset(gsea_kegg_res@result, NES > 0 & p.adjust < 0.05)
qs::qsave(gsea_kegg_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/03_gsea_kegg.qs")

gsea_kegg_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/03_gsea_kegg.qs")
prompt <- interpret_tool(gsea_kegg_res, database = "KEGG", diff_genes, pathway_num = 50, contact_LLM = FALSE,
                         context = "786-0 cells were treated with either a control DMSO vehicle or 10um Everolimus, an mTOR inhibitor, for 4 hours.")
writeLines(prompt, "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/prompt/gsea_kegg.txt")

##### plot #####
gsea_kegg_dat <- as.data.frame(gsea_kegg_res) %>% dplyr::arrange(desc(NES)) %>% dplyr::slice_head(n = 50)
p <- KEGGStratifyPlot(enrich_data = gsea_kegg_dat,
                      kegg_df = kegg_df,
                      enrich_score = "NES",
                      pval = "p.adjust")
print(p)
ggsave("analysis/figure/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/kegg_classified.pdf", p, width = 11, height = 7)

#### WP ####
ora_wp_res <- enrichWP(gene         = entrez, 
                       organism     = "Homo sapiens",
                       pvalueCutoff = 0.05,
                       universe     = univers_entrez,
                       minGSSize    = 0,
                       maxGSSize    = 1000000) 
ora_wp_res <- setReadable(ora_wp_res, 
                          OrgDb = "org.Hs.eg.db",
                          keyType = "ENTREZID")
ora_wp_res@result <- subset(ora_wp_res@result, p.adjust < 0.05)
ora_wp_dat <- as.data.frame(ora_wp_res) %>% dplyr::arrange(desc(FoldEnrichment))
qs::qsave(ora_wp_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/04_ora_wp.qs")


gsea_wp_res <- gseWP(geneList     = entrez_list, 
                     organism     = "Homo sapiens",
                     minGSSize    = 0,
                     maxGSSize    = 1000000,
                     pvalueCutoff = 0.05,
                     verbose      = TRUE)
gsea_wp_res <- setReadable(gsea_wp_res,
                           OrgDb = "org.Hs.eg.db",
                           keyType = "ENTREZID")
gsea_wp_res@result <- subset(gsea_wp_res@result, NES > 0 & p.adjust < 0.05)
gsea_wp_dat <- as.data.frame(gsea_wp_res) %>% arrange(desc(NES))
qs::qsave(gsea_wp_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/04_gsea_wp.qs")

gsea_wp_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/04_gsea_wp.qs")
prompt <- interpret_tool(gsea_wp_res, database = "WikiPathways", diff_genes, pathway_num = NULL, contact_LLM = FALSE,
                         context = "786-0 cells were treated with either a control DMSO vehicle or 10um Everolimus, an mTOR inhibitor, for 4 hours.")
writeLines(prompt, "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/prompt/gsea_wp.txt")

#### Reactome ####
ora_reactome_res <- enrichPathway(gene = entrez, 
                                  organism = "human",
                                  pvalueCutoff = 0.05,
                                  universe = univers_entrez,
                                  minGSSize    = 0,
                                  maxGSSize    = 1000000)
ora_reactome_res <- setReadable(ora_reactome_res,
                                OrgDb = "org.Hs.eg.db",
                                keyType = "ENTREZID")
ora_reactome_res@result <- subset(ora_reactome_res@result, p.adjust < 0.05)
ora_reactome_dat <- as.data.frame(ora_reactome_res) %>% dplyr::arrange(desc(FoldEnrichment))
qs::qsave(ora_reactome_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/05_ora_reactome.qs")

gsea_reactome_res<- gsePathway(geneList = entrez_list, 
                               organism = "human",
                               minGSSize    = 0,
                               maxGSSize    = 1000000,
                               pvalueCutoff = 0.05,
                               verbose = TRUE)
gsea_reactome_res <- setReadable(gsea_reactome_res,
                                 OrgDb = "org.Hs.eg.db",
                                 keyType = "ENTREZID")
gsea_reactome_res@result <- subset(gsea_reactome_res@result, NES > 0 & p.adjust < 0.05)
reactome_res_dat <- as.data.frame(gsea_reactome_res) %>% arrange(desc(NES))
qs::qsave(gsea_reactome_res, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/05_gsea_reactome.qs")

gsea_reactome_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/05_gsea_reactome.qs")
prompt <- interpret_tool(gsea_reactome_res, database = "Reactome", diff_genes, pathway_num = 53, contact_LLM = FALSE,
                         context = "786-0 cells were treated with either a control DMSO vehicle or 10um Everolimus, an mTOR inhibitor, for 4 hours.")
writeLines(prompt, "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/prompt/gsea_reactome.txt")

#### Determine ERGs ####
ERGs <- deg_res %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange > 1) %>% 
  dplyr::filter(log2(baseMean) > 10) %>% 
  pull(gene)
qs::qsave(ERGs, file = "analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/06_ERGs.qs")

data <- deg_res
logFCfilter = 1
data$is_ERG <- ifelse(data$gene %in% ERGs, "ERGs", "Non-ERGs")
y_limit <- max(abs(data$log2FoldChange), na.rm = TRUE)*1.1

p <- ggplot(data = data, 
            aes(x = log2(baseMean), 
                y = log2FoldChange, 
                color = is_ERG)) +  # 修改颜色映射为is_ERG列
  geom_point(size  = 0.6)+
  scale_color_manual(values = c("red4", "grey80"))+  # ERG基因红色，其他黑色
  scale_y_continuous(limits = c(-y_limit, y_limit))+
  labs(y = "log2 (Fold Change)", 
       x = "log2 (Base Mean)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = c(logFCfilter,-logFCfilter),
             lty = 2,
             lwd = 1)+
  geom_hline(yintercept = 0,
             lwd = 1.2)+
  theme_bw()+
  theme(panel.border     = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black")) +
  theme(legend.position = "right")+  # 显示图例以便区分
  labs(color = "Gene Type") +  # 修改图例标题
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        aspect.ratio = 1)

print(p)
ggsave(filename = "analysis/figure/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/DESeq2_Maplot_ERGs.pdf", p, height = 3.96, width = 4.8)
