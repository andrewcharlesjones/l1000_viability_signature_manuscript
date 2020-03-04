library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)
library(reshape2)
library(GSEABase)


plot_dir <- "~/Documents/cancer_data_science/l1000_via_sig_ms/figures/plots/"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)

l1000_shrna_data <- l1000.analysis::load_l1000_shrna_data()


### COMPUTE PERT-SPECIFIC SIGS --------

cpd_viability_sigs <- l1000.analysis::compute_viability_signature_all_cpds(data_list = l1000_data, sens_column = "auc_avg")

shrna_via_sigs <- l1000.analysis::compute_viability_signature_all_cpds(data_list = l1000_shrna_data, sens_column = "shrna_viability")


### GET MATRICES OF VIA SIGS ------


cpd_viability_sigs_melted <- melt(cpd_viability_sigs)
cpd_viability_related_sigs <- cpd_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "logFC")
cpd_viability_related_sigs_cast <- reshape2::acast(cpd_viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()


shrna_viability_sigs_melted <- melt(shrna_via_sigs)
shrna_viability_related_sigs <- shrna_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "logFC")
shrna_viability_related_sigs_cast <- reshape2::acast(shrna_viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()


stopifnot(all.equal(rownames(cpd_viability_related_sigs_cast), rownames(shrna_viability_related_sigs_cast)))


combined_sigs <- cbind(cpd_viability_related_sigs_cast, shrna_viability_related_sigs_cast)


sig_cors <- cor(combined_sigs)
diag(sig_cors) <- NA

pheatmap::pheatmap(sig_cors)


### HEATMAP SMALLER SUBSET OF PERTS TO SHOW CONTRAST BETWEEN CLUSTERS ------

high_cor_idx <- which(apply(sig_cors, 1, max, na.rm=T) > 0.5)
sig_cors_extreme <- sig_cors[high_cor_idx, high_cor_idx]
pheatmap::pheatmap(sig_cors_extreme, height = 8, width = 8, filename = file.path(plot_dir, "via_sig_heatmap_small_subset.png"))



### TSNE OF VIA SIGS ------

combined_sigs_t <- combined_sigs %>% t()
combined_sigs_t_nodup <- combined_sigs_t[!duplicated(rownames(combined_sigs_t)),]


tsne_out <- Rtsne::Rtsne(combined_sigs_t_nodup[,lm_genes])

tsne_coords <- tsne_out$Y %>% 
  as.data.frame() %>% 
  set_colnames(c("t1", "t2")) %>% 
  set_rownames(rownames(combined_sigs_t_nodup))

tsne_coords %>% ggplot(aes(t1, t2)) +
  geom_point() + 
  cdsr::theme_Publication()

repurp_data <- l1000.analysis::load_repurposing_cpd_annotations()

tsne_coords_annot <- merge(tsne_coords, repurp_data[,c("Name", "MOA")], by.x = "row.names", by.y = "Name", all.x = T)

tsne_coords_annot %>% ggplot(aes(t1, t2, color = tsne_coords_annot$MOA == "BCL inhibitor")) +
  geom_point() + 
  cdsr::theme_Publication()


