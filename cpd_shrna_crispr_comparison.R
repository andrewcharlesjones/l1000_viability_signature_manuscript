library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)
library(reshape2)

plot_dir <- "/Users/andrewjones/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots"


### LOAD DATA ----------

# compounds
l1000_cpd_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data_cpd <- cbind(l1000_cpd_data$expression, l1000_cpd_data$sensitivity, l1000_cpd_data$metadata)

merged_data_cpd_collapsed <- merged_data_cpd %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)

# shRNA
l1000_shrna_data <- l1000.analysis::load_l1000_shrna_data()
merged_data_shrna <- cbind(l1000_shrna_data$expression, l1000_shrna_data$sensitivity, l1000_shrna_data$metadata)

merged_data_shrna_collapsed <- merged_data_shrna %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)

# CRISPR
l1000_crispr_data <- l1000.analysis::load_l1000_crispr_data()
merged_data_crispr <- cbind(l1000_crispr_data$expression, l1000_crispr_data$sensitivity, l1000_crispr_data$metadata)

merged_data_crispr_collapsed <- merged_data_crispr %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)



### COMPUTE GLOBAL SIGNATURES -------

# cpds
cpd_global_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_cpd_data, global = T, sens_column = "auc_avg")

# shRNA
shrna_global_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_shrna_data, global = T, sens_column = "shrna_viability")

# CRISPR
crispr_global_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_crispr_data, global = T, sens_column = "crispr_viability")



### HEATMAP SIGNATURES --------

relevant_cols <- c("logFC", "Gene")
cpd_and_sh_sigs <- merge(cpd_global_via_sig$viability_related[,relevant_cols], shrna_global_via_sig$viability_related[,relevant_cols], by = "Gene",
                         suffixes = c("_cpd", "_shrna"))
all_three_sigs <- merge(cpd_and_sh_sigs, crispr_global_via_sig$viability_related[,relevant_cols] %>% 
                          dplyr::rename(logFC_crispr = logFC),
                        by = "Gene") %>% 
  column_to_rownames("Gene")

num_top_genes <- 30
top_genes <- rownames(all_three_sigs)[order(-abs(rowMeans(all_three_sigs)))][1:num_top_genes]


pheatmap::pheatmap(all_three_sigs[top_genes,] %>% scale() %>% t(), 
                   labels_row = c("RNAi", "Compounds", "CRISPR"), 
                   cellheight = 25,
                   cellwidth = 10,
                   filename = file.path(plot_dir, "cpd_shrna_via_sig_comparison.png"))

pheatmap::pheatmap(all_three_sigs[top_genes,] %>% scale() %>% t(), 
                   labels_row = c("RNAi", "Compounds", "CRISPR"), height = 5, width = 15,
                   filename = file.path(plot_dir, "cpd_shrna_via_sig_comparison.png"), fontsize = 20)

cor(all_three_sigs)

gene_info <- l1000.analysis::load_gene_info()
lm_genes <- gene_info$pr_gene_symbol[gene_info$pr_is_lm == 1]

cor(all_three_sigs[lm_genes,], use = "pairwise.complete.obs")



### LOOK AT RELATIONSHIP BETWEEN LM AND INFERRED GENES --------

## logFC
relevant_cols <- c("logFC", "Gene")

cpd_and_sh_sigs <- merge(cpd_global_via_sig$viability_related[,relevant_cols], shrna_global_via_sig$viability_related[,relevant_cols], by = "Gene",
                         suffixes = c("_cpd", "_shrna"))
all_three_sigs <- merge(cpd_and_sh_sigs, crispr_global_via_sig$viability_related[,relevant_cols] %>% 
                          dplyr::rename(logFC_crispr = logFC),
                        by = "Gene")

all_three_sigs_melted <- melt(all_three_sigs)

all_three_sigs_melted <- merge(all_three_sigs_melted, gene_info[,c("pr_gene_symbol", "pr_is_lm")], by.x = "Gene", by.y = "pr_gene_symbol", all.x = T)

p1 <- all_three_sigs_melted %>% ggplot(aes(variable %>% as.factor(), value %>% abs(), color = pr_is_lm %>% as.factor())) + 
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  labs(x = "", y = "logFC") + 
  l1000.analysis::theme_scatterplot() + 
  guides(color=guide_legend(title="Landmark gene"))


## t-statistics
# relevant_cols <- c("t", "Gene")
# 
# cpd_and_sh_sigs <- merge(cpd_global_via_sig$viability_related[,relevant_cols], shrna_global_via_sig$viability_related[,relevant_cols], by = "Gene",
#                          suffixes = c("_cpd", "_shrna"))
# all_three_sigs <- merge(cpd_and_sh_sigs, crispr_global_via_sig$viability_related[,relevant_cols] %>% 
#                           dplyr::rename(t_crispr = t),
#                         by = "Gene")
# 
# all_three_sigs_melted <- melt(all_three_sigs)
# 
# all_three_sigs_melted <- merge(all_three_sigs_melted, gene_info[,c("pr_gene_symbol", "pr_is_lm")], by.x = "Gene", by.y = "pr_gene_symbol", all.x = T)
# 
# p2 <- all_three_sigs_melted %>% ggplot(aes(variable %>% as.factor(), value %>% abs(), color = pr_is_lm %>% as.factor())) + 
#   geom_boxplot() + 
#   cdsr::theme_Publication() + 
#   labs(x = "", y = "t-stat") + 
#   l1000.analysis::theme_scatterplot() + 
#   guides(color=guide_legend(title="Landmark gene"))
# 
# 
# p_grid <- cowplot::plot_grid(p1, p2)
# ggsave(file.path(plot_dir, "lm_vs_inferred_via_sig_strength.png"), p_grid, height = 10, width = 20)













