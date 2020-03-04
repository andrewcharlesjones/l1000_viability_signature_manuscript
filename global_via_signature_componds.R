library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)
library(ggrepel)

plot_dir <- "/Users/andrewjones/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data <- cbind(l1000_data$expression, l1000_data$sensitivity, l1000_data$metadata)

merged_data_collapsed <- merged_data %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)


### PLOT EXAMPLE EXPRESSION-VIABILITY RELATIONSHIPS -----------

plot_global_via_relationship <- function(gene_name) {
  
  merged_data_collapsed %>% ggplot(aes_string("auc_avg", gene_name)) + 
    geom_point() + 
    geom_smooth(method = "lm", size = 4) + 
    # labs(x = "Sensitivity", y = latex2exp::TeX(paste('$\\Delta$', gene_name, sep = " "))) + 
    labs(x = "Sensitivity", y = paste(gene_name, "DE", sep = " ")) + 
    cdsr::theme_Publication() + 
    l1000.analysis::theme_scatterplot() + 
    theme(axis.title = element_text(face = "plain"))
  
}


p1 <- plot_global_via_relationship("CDC20")
p2 <- plot_global_via_relationship("GADD45A")

p_grid <- cowplot::plot_grid(p1, p2)
ggsave(file.path(plot_dir, "gene_sensitivity_relationship_examples.png"), p_grid, height = 5, width = 10)






### COMPUTE GLOBAL SIG ----------

global_viability_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_avg")


apoptosis_genes <- l1000.analysis::get_gene_set_genes(gsc = l1000.analysis::load_gsc_data(), collection_name = "hallmark", gene_set_name = "HALLMARK_APOPTOSIS")
cc_genes <- c(
  l1000.analysis::get_gene_set_genes(gsc = l1000.analysis::load_gsc_data(), collection_name = "hallmark", gene_set_name = "HALLMARK_E2F_TARGETS"),
  l1000.analysis::get_gene_set_genes(gsc = l1000.analysis::load_gsc_data(), collection_name = "hallmark", gene_set_name = "HALLMARK_G2M_CHECKPOINT")
)

df_for_plotting <- global_viability_sig$viability_related
gene_class_list <- list()

gene_class_list <- llply(df_for_plotting$Gene, function(x) {
  if (x %in% apoptosis_genes) {
    return("Apoptosis")
  } else if (x %in% cc_genes) {
    return("Cell cycle")
  } else {
    return(NA)
  }
}) %>% as.character()
df_for_plotting['gene_class'] <- gene_class_list


p <- df_for_plotting %>% 
  ggplot(aes(logFC, -log10(P.Value), color = gene_class)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  scale_color_manual(labels = c("Apoptosis", "Cell cycle", ""), values=c("red", "blue", "gray")) + 
  geom_label_repel(data = rbind(
    df_for_plotting[df_for_plotting$gene_class == "Apoptosis",] %>% arrange(-logFC) %>% head(10),
    df_for_plotting[df_for_plotting$gene_class == "Cell cycle",] %>% arrange(logFC) %>% head(10)),
                   aes(label = Gene),
    size = 10) + 
  labs(x = "Viability signature coefficient", y = "-log10(p-value)") + 
  theme(legend.title=element_blank())

ggsave(file.path(plot_dir, "global_via_sig_volcano.png"), p, height = 10, width = 10)





                     


### GSEA ON GLOBAL SIG -------

gsc <- l1000.analysis::load_gsc_data()

gsea_hallmark <- cdsr::run_fGSEA(gene_stat = global_viability_sig$viability_related$logFC %>% set_names(global_viability_sig$viability_related$Gene), 
                gsc = gsc[["hallmark"]], 
                perm_type = "gene")

gsea_hallmark %>% arrange(NES) %>% head()
gsea_hallmark %>% arrange(-NES) %>% head()



# down table
gsea_table_df_down <- gsea_hallmark[,c("NES", "pval", "pathway")] %>% arrange(NES) %>% head(5) %>% arrange(pval)
gsea_table_df_down <- gsea_table_df_down[,c("pathway", "pval")]
gsea_table_df_down$pathway <- stringr::str_replace_all(gsea_table_df_down$pathway, pattern = "HALLMARK_", replacement = "")
gsea_table_df_down$pathway <- stringr::str_replace_all(gsea_table_df_down$pathway, pattern = "_", replacement = " ")
gsea_table_df_down$pval <- format(round(gsea_table_df_down$pval, 5), scientific = T)

dev.off()
png(file.path(plot_dir, "gsea_down.png"), height = 100*nrow(gsea_table_df_down), width = 400*ncol(gsea_table_df_down))

gridExtra::grid.table(gsea_table_df_down)
dev.off()


# up table
gsea_table_df_up <- gsea_hallmark[,c("NES", "pval", "pathway")] %>% arrange(-NES) %>% head(5) %>% arrange(pval)
gsea_table_df_up <- gsea_table_df_up[,c("pathway", "pval")]
gsea_table_df_up$pathway <- stringr::str_replace_all(gsea_table_df_up$pathway, pattern = "HALLMARK_", replacement = "")
gsea_table_df_up$pathway <- stringr::str_replace_all(gsea_table_df_up$pathway, pattern = "_", replacement = " ")
gsea_table_df_up$pval <- format(round(gsea_table_df_up$pval, 5), scientific = T)

dev.off()
png(file.path(plot_dir, "gsea_up.png"), height = 100*nrow(gsea_table_df_up), width = 400*ncol(gsea_table_df_up))

gridExtra::grid.table(gsea_table_df_up)
dev.off()











### COMPARE GLOBAL VIA SIG OF DRUG CATEGORIES -------

# targeted cancer
targeted_cancer_idx <- l1000_data$metadata$pert_iname %in% corsello.drug.categories$name[corsello.drug.categories$drug_category == "targeted cancer"]
l1000_data_targeted_cancer <- list(expression = l1000_data$expression[targeted_cancer_idx,],
                                   sensitivity = l1000_data$sensitivity[targeted_cancer_idx,],
                                   metadata = l1000_data$metadata[targeted_cancer_idx,])

targeted_cancer_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data_targeted_cancer, global = T, sens_column = "auc_avg")


# chemo
targeted_cancer_idx <- l1000_data$metadata$pert_iname %in% corsello.drug.categories$name[corsello.drug.categories$drug_category == "chemo"]
l1000_data_targeted_cancer <- list(expression = l1000_data$expression[targeted_cancer_idx,],
                                   sensitivity = l1000_data$sensitivity[targeted_cancer_idx,],
                                   metadata = l1000_data$metadata[targeted_cancer_idx,])

chemo_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data_targeted_cancer, global = T, sens_column = "auc_avg")

# noncancer
targeted_cancer_idx <- l1000_data$metadata$pert_iname %in% corsello.drug.categories$name[corsello.drug.categories$drug_category == "noncancer"]
l1000_data_targeted_cancer <- list(expression = l1000_data$expression[targeted_cancer_idx,],
                                   sensitivity = l1000_data$sensitivity[targeted_cancer_idx,],
                                   metadata = l1000_data$metadata[targeted_cancer_idx,])

noncancer_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data_targeted_cancer, global = T, sens_column = "auc_avg")


# look at diffs

gene_names <- noncancer_via_sig$viability_related$Gene

three_sigs <- data.frame(targeted_cancer = (targeted_cancer_via_sig$viability_related %>% column_to_rownames("Gene"))[gene_names,]$logFC,
                         chemo = (chemo_via_sig$viability_related %>% column_to_rownames("Gene"))[gene_names,]$logFC,
                         noncancer = (noncancer_via_sig$viability_related %>% column_to_rownames("Gene"))[gene_names,]$logFC,
                         row.names =gene_names
)

three_sigs_cors <- cor(three_sigs)
diag(three_sigs_cors) <- NA
pheatmap::pheatmap(three_sigs_cors)

three_sigs_df <- three_sigs %>% rownames_to_column("Gene")

gene_type <- llply(three_sigs_df$Gene, function(x) {
  if (x %in% apoptosis_genes) {
    return("apoptosis") 
  } else if (x %in% cell_cycle_genes) {
    return("cell cycle") 
  } else {
    return("other")
  }
}) %>% as.character()
three_sigs_df$gene_type <- gene_type

three_sigs_df %>% ggplot(aes(targeted_cancer, chemo, color = gene_type)) + 
  geom_point() + 
  geom_text_repel(data = unique(rbind(three_sigs_df %>% arrange(-abs(chemo)) %>% head(20),
                                      three_sigs_df %>% arrange(-abs(targeted_cancer)) %>% head(20))),
                  aes(label = Gene)) + 
  cdsr::theme_Publication() + 
  scale_color_manual(values=c("red", "blue", "gray"))


chemo_gsea <- cdsr::run_fGSEA(gsc = gsc[["hallmark"]], gene_stat = three_sigs_df$chemo %>% set_names(three_sigs_df$Gene), perm_type = "gene")
chemo_gsea %>% arrange(-abs(NES)) %>% head(10)

targeted_gsea <- cdsr::run_fGSEA(gsc = gsc[["hallmark"]], gene_stat = three_sigs_df$targeted_cancer %>% set_names(three_sigs_df$Gene), perm_type = "gene")
targeted_gsea %>% arrange(-abs(NES)) %>% head(10)
















### COMPARE GLOBAL VIA SIG OF VIABILITY DATASETS -------

global_viability_sig_prism <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_prism")
global_viability_sig_gdsc <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_gdsc")
global_viability_sig_ctrp <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_ctrp")

prism_sig <- global_viability_sig_prism$viability_related %>% set_rownames(global_viability_sig_prism$viability_related$Gene)
gdsc_sig <- global_viability_sig_gdsc$viability_related %>% set_rownames(global_viability_sig_gdsc$viability_related$Gene)
ctrp_sig <- global_viability_sig_ctrp$viability_related %>% set_rownames(global_viability_sig_ctrp$viability_related$Gene)

all_global_sigs_df <- data.frame(
  PRISM = prism_sig[lm_genes, "logFC"],
  GDSC = gdsc_sig[lm_genes, "logFC"],
  CTRP = ctrp_sig[lm_genes, "logFC"]
)

all_global_sig_cors <- cor(all_global_sigs_df)
diag(all_global_sig_cors) <- NA
pheatmap::pheatmap(all_global_sig_cors, 
                   cluster_cols = F, 
                   cluster_rows = F, 
                   breaks=seq(0, 1, by = .01), 
                   display_numbers = T, fontsize = 15, filename = "~/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots/global_sig_comparison_across_datasets.png")




### COMPARE R^2 BETWEEN LANDMARK AND NON-LANDMARK GENES -----

rsq <- function (x, y) cor(x, y, use = "pairwise.complete.obs") ^ 2

cpd_data_with_meta <- cbind(l1000_data$expression,
                            l1000_data$sensitivity,
                            l1000_data$metadata)

cpd_data_collapsed <- cpd_data_with_meta %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)


l1k_gene_info <- l1000.analysis::load_gene_info()

rsq_list_lm <- list()
for (curr_gene in lm_genes) {
  curr_rsq <- rsq(cpd_data_collapsed[,curr_gene], cpd_data_collapsed$auc_avg)
  rsq_list_lm <- c(rsq_list_lm, curr_rsq)
}

inferred_genes <- l1k_gene_info$pr_gene_symbol[l1k_gene_info$pr_is_lm == 0]
rsq_list_inf <- list()
for (curr_gene in inferred_genes) {
  curr_rsq <- rsq(cpd_data_collapsed[,curr_gene], cpd_data_collapsed$auc_avg)
  rsq_list_inf <- c(rsq_list_inf, curr_rsq)
}

r2_df <- data.frame(
  r2 = c(rsq_list_lm %>% as.double(), rsq_list_inf %>% as.double()),
  gene_type = c(rep("Landmark", length(rsq_list_lm)), rep("Inferred", length(rsq_list_inf)))
)

r2_df %>% ggplot(aes(gene_type, r2)) + 
  geom_boxplot()


