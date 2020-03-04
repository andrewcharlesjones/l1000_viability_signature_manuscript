library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)
library(reshape2)
library(GSEABase)
library(ggrepel)

plot_dir <- "~/Documents/cancer_data_science/l1000_via_sig_ms/figures/plots/"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data <- cbind(l1000_data$expression, l1000_data$sensitivity, l1000_data$metadata)



### COMPUTE GLOBAL SIG FOR EACH TIMEPOINT ----------

tp24_idx <- l1000_data$metadata$pert_time == 24
tp24_data_list <- list(
  expression = l1000_data$expression[tp24_idx,],
  sensitivity = l1000_data$sensitivity[tp24_idx,],
  metadata = l1000_data$metadata[tp24_idx,]
)

tp6_idx <- l1000_data$metadata$pert_time == 6
tp6_data_list <- list(
  expression = l1000_data$expression[tp6_idx,],
  sensitivity = l1000_data$sensitivity[tp6_idx,],
  metadata = l1000_data$metadata[tp6_idx,]
)

global_viability_sig_24hr <- l1000.analysis::compute_viability_signature(data_list = tp24_data_list, global = T, sens_column = "auc_avg")
global_viability_sig_6hr <- l1000.analysis::compute_viability_signature(data_list = tp6_data_list, global = T, sens_column = "auc_avg")



both_tp_sigs <- merge(global_viability_sig_24hr$viability_related, global_viability_sig_6hr$viability_related, by = "Gene", suffixes = c("_24hr", "_6hr"))

p <- both_tp_sigs %>% ggplot(aes(logFC_6hr, logFC_24hr)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_abline() + 
  labs(x = "6-hour signature", y = "24-hour signature") + 
  geom_label_repel(data = unique(rbind(both_tp_sigs %>% arrange(-abs(logFC_24hr)) %>% head(6),
                                       both_tp_sigs %>% arrange(-abs(logFC_6hr)) %>% head(6))),
                   aes(label = Gene),
                   size = 10)
ggsave(file.path(plot_dir, "timepoint_comparison.png"), p, height = 10, width = 11)











### LOOK FOR DIFFERENCES BETWEEN THE TIMEPOINTS -------

tp_sigs_diffs <- (both_tp_sigs$logFC_24hr %>% scale()) - (both_tp_sigs$logFC_6hr %>% scale())

both_tp_sigs['diff'] <- tp_sigs_diffs


both_tp_sigs_scaled <- both_tp_sigs %>% 
  mutate(logFC_6hr = logFC_6hr %>% scale(),
         logFC_24hr = logFC_24hr %>% scale())

both_tp_sigs_scaled %>% ggplot(aes(logFC_6hr, logFC_24hr, color = diff %>% abs())) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_abline() + 
  labs(x = "6-hour signature", y = "24-hour signature") + 
  geom_label_repel(data = unique(rbind(both_tp_sigs_scaled %>% arrange(-abs(diff)) %>% head(20))),
                                       # both_tp_sigs %>% arrange(-abs(logFC_6hr)) %>% head(10))),
                   aes(label = Gene))


gsc_data <- l1000.analysis::load_gsc_data()
gsea_out <- cdsr::run_fGSEA(gene_stat = both_tp_sigs_scaled$diff %>% set_names(both_tp_sigs_scaled$Gene), gsc = gsc_data[["GO_biological_process"]], perm_type = "gene")
gsea_out %>% arrange(pval) %>% head(15)


