library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)

plot_dir <- "~/Documents/cancer_data_science/l1000_via_sig_ms/figures/plots/"


### LOAD DATA --------

l1000_shrna_data <- l1000.analysis::load_l1000_shrna_data()
full_shrna_data <- cbind(l1000_shrna_data$expression, l1000_shrna_data$metadata, l1000_shrna_data$sensitivity)

full_shrna_data_collapsed <- full_shrna_data %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)


### GLOBAL SH GLOBAL VIABILITY SIG ---------

shrna_global_via_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_shrna_data, global = T, sens_column = "shrna_viability")

### COMPUTE SH VIABILITY SIGS ---------

shrna_viability_sigs <- l1000.analysis::compute_viability_signature_all_cpds(data_list = l1000_shrna_data, sens_column = "shrna_viability")


shrna_viability_sigs_melted <- melt(shrna_viability_sigs)
shrna_viability_related_sigs <- shrna_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "t")
shrna_viability_related_sigs_cast <- reshape2::acast(shrna_viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()

via_sig_cors <- cor(shrna_viability_related_sigs_cast)
diag(via_sig_cors) <- NA
pheatmap::pheatmap(via_sig_cors, fontsize = 5)


### PLOT VOLCANO OF SIGS --------

shrna_to_plot <- "MAP2K1"
via_sig <- shrna_viability_sigs[[shrna_to_plot]]$viability_related

p <- via_sig %>% ggplot(aes(logFC, -log10(P.Value))) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  labs(x = "Viability-related signature") + 
  geom_label_repel(data = via_sig %>% arrange(P.Value) %>% head(15),
                   aes(label = Gene),
                   size = 10) + 
  ggtitle(paste(shrna_to_plot, "KD", sep = " "))
p
ggsave(file.path(plot_dir, "map2k1_via_sig_volcano.png"), p, height = 10, width = 10)




p <- l1000.analysis::make_via_sig_volcano(via_sig_output = shrna_viability_sigs[[shrna_to_plot]], cpd_name = shrna_to_plot)
p




### COMPUTE CORRELATIONS BETWEEN RELATED AND UNRELATED SIGS --------

sig_correlations <- ldply(shrna_viability_sigs, function(x) {
  combined_data <- merge(x$viability_related, x$viability_unrelated,
                         by = "Gene",
                         suffixes = c("_slope", "_intercept"))
  return(cor(combined_data$logFC_slope, combined_data$logFC_intercept))
})
sig_correlations %>% arrange(abs(V1)) %>% head(15)


### PLOT VIABILITY-RELATED AND -UNRELATED AGAINST EACH OTHER --------

l1000.analysis::plot_via_sig_comparison(via_sig_output = shrna_viability_sigs$MAP2K1)


### PLOT HEATMAP OF EXPRESSION RESPONSE ------

gene_to_heatmap <- "MAP2K1"
l1000.analysis::plot_pert_data_heatmap(data_list = l1000_shrna_data, via_sig_list = shrna_viability_sigs, cpd_name = gene_to_heatmap, 
                                       sens_columns = c("shrna_viability"), sort_column = "shrna_viability", cellheight = 10, cellwidth = 10, remove_via_nas = T,
                                       filename = file.path(plot_dir, "map2k1_heatmap.png"))


### COMPARE CPD AND SHRNA SIGS ---------

shared_pathway_genes <- c("MYC", "DUSP6", "FOSL1", "IER3")

p <- l1000.analysis::plot_cpd_and_shrna_sigs(cpd_via_sig_df = cpd_viability_sigs, shrna_via_sig_df = shrna_viability_sigs, cpd_name = "selumetinib", shrna_name = "MAP2K1", shared_pathway_genes = shared_pathway_genes) + 
  scale_color_manual(labels = c("", "MAPK pathway"), values = c("gray", "red"))
p
ggsave(file.path(plot_dir, "cpd_shrna_via_sig_comparison.png"), p, height = 10, width = 10)



