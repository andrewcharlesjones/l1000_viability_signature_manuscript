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



plot_dir <- "~/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots/"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data <- cbind(l1000_data$expression, l1000_data$sensitivity, l1000_data$metadata)

gene_info <- l1000.analysis::load_gene_info()
lm_genes <- gene_info$pr_gene_symbol[gene_info$pr_is_lm == 1]


### COMPUTE VIABILITY SIGNATURES FOR EACH PERT --------

cpd_viability_sigs <- l1000.analysis::compute_viability_signature_all_cpds_control_dose_and_time(data_list = l1000_data, sens_column = "auc_avg")

cpd_viability_sigs_melted <- melt(cpd_viability_sigs)
viability_related_sigs <- cpd_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "t")
viability_related_sigs_cast <- reshape2::acast(viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()


l1000_shrna_data <- l1000.analysis::load_l1000_shrna_data()

shrna_viability_sigs <- l1000.analysis::compute_viability_signature_all_cpds_control_dose_and_time(data_list = l1000_shrna_data, sens_column = "shrna_viability")


shrna_viability_sigs_melted <- melt(shrna_viability_sigs)
shrna_viability_related_sigs <- shrna_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "t")
shrna_viability_related_sigs_cast <- reshape2::acast(shrna_viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()


all_viability_sigs <- merge(viability_related_sigs_cast, shrna_viability_related_sigs_cast, by = "row.names")
all_viability_sigs %<>% column_to_rownames("Row.names")



annot_df <- ldply(colnames(all_viability_sigs), function(x) {
  has_sanger <- !is.na(merged_data[merged_data$pert_iname == x, "auc_gdsc"][1])
  has_prism <- !is.na(merged_data[merged_data$pert_iname == x, "auc_prism"][1])
  has_ctrp <- !is.na(merged_data[merged_data$pert_iname == x, "auc_ctrp"][1])
  
  return(data.frame(
    has_sanger = has_sanger %>% as.integer(),
    has_prism = has_prism %>% as.integer(),
    has_ctrp = has_ctrp %>% as.integer()
  ))
}) %>% 
  set_rownames(colnames(all_viability_sigs))


via_sig_cors <- cor(all_viability_sigs[lm_genes,])
# diag(via_sig_cors) <- NA
pheatmap::pheatmap(via_sig_cors, annotation_col = annot_df, annotation_row = annot_df)
pheatmap::pheatmap(via_sig_cors, filename = "~/Desktop/a.png", height = 10, width = 10, fontsize = 3)


### DEFINE TWO BIG CLUSTERS -------

clust_out <- cutree(hclust(dist(via_sig_cors)), 12)
weird_cluster_perts <- names(clust_out)[clust_out == clust_out[names(clust_out) == "BAX-channel-blocker"]]
clust_out2 <- cutree(hclust(dist(via_sig_cors)), 4)
main_cluster_perts <- names(clust_out2)[clust_out2 == clust_out2[names(clust_out2) == "nutlin-3"]]

both_cluster_sigs <- data.frame(
  weird_cluster = all_viability_sigs[,weird_cluster_perts] %>% rowMeans(),
  main_cluster = all_viability_sigs[,main_cluster_perts] %>% rowMeans()
)
both_cluster_sigs %<>% mutate(Gene = rownames(both_cluster_sigs))

apoptosis_genes <- read.table("~/Documents/cancerdatascience/l1000_via_sig_ms/apoptosis_genes.txt") %>% 
  .$V2 %>% as.character()
apoptosis_genes <- apoptosis_genes[apoptosis_genes != "GeneName"]

cell_cycle_genes <- read.table("~/Documents/cancerdatascience/l1000_via_sig_ms/cell_cycle_genes.txt") %>% 
  .$V2 %>% as.character()
cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes != "GeneName"]

proliferation_genes <- read.table("~/Documents/cancerdatascience/l1000_via_sig_ms/proliferation_genes.txt") %>% 
  .$V2 %>% as.character()
proliferation_genes <- proliferation_genes[proliferation_genes != "GeneName"]

both_cluster_sigs$gene_class <- llply(both_cluster_sigs$Gene, function(x) {
  if (x %in% apoptosis_genes) {
    return("apoptosis")
  } else if (x %in% cell_cycle_genes) {
    return("cell cycle")
  } else if (x %in% proliferation_genes) {
    return("proliferation")
  } else {
    return("")
  }
}) %>% as.character()

p <- both_cluster_sigs %>% ggplot(aes(weird_cluster, main_cluster, color = gene_class)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  scale_color_manual(values=c("gray", "red", "blue", "green"))
ggExtra::ggMarginal(p, groupColour = T)










### HEATMAP VIABILITY SIGS -------

top_genes <- c(
  rownames(all_viability_sigs)[order(all_viability_sigs[,weird_cluster_perts] %>% rowMeans())[1:10]],
  rownames(all_viability_sigs)[order(-all_viability_sigs[,weird_cluster_perts] %>% rowMeans())[1:10]]
)

annot_df <- ldply(colnames(all_viability_sigs), function(x) {
  if (x %in% weird_cluster_perts) {
    return("Cluster 2")
  } else if (x %in% main_cluster_perts) {
    return("Cluster 1")
  } else {
    return("")
  }
}) %>% 
  set_colnames("cluster") %>% 
  set_rownames(colnames(all_viability_sigs))

sorted_idx <- order(annot_df)

pheatmap::pheatmap(all_viability_sigs[top_genes,sorted_idx], annotation_col = annot_df[sorted_idx,,drop = F], cluster_cols = F)
pheatmap::pheatmap(all_viability_sigs[top_genes,sorted_idx], annotation_col = annot_df[sorted_idx,,drop = F], cluster_cols = F, show_colnames = F,
                   filename = file.path(plot_dir, "viability_sig_heatmap_clustering.png"), cellheight = 10)



### PLOT ONE-BY-ONE EXAMPLES FROM EACH CLUSTER ------

via_sigs_plotting <- all_viability_sigs %>% mutate(Gene = rownames(all_viability_sigs))
via_sigs_plotting$gene_class <- llply(via_sigs_plotting$Gene, function(x) {
  if (x %in% apoptosis_genes) {
    return("apoptosis")
  } else if (x %in% cell_cycle_genes) {
    return("cell cycle")
  } else if (x %in% proliferation_genes) {
    return("proliferation")
  } else {
    return("")
  }
}) %>% as.character()

p <- via_sigs_plotting[lm_genes,] %>% ggplot(aes(`YM-155`, `BAX-channel-blocker`, color = gene_class)) + 
  geom_point() + 
  cdsr::theme_Publication()
ggExtra::ggMarginal(p, groupColour = T)


e2f_genes <- l1000.analysis::get_gene_set_genes(gsc = gsc_data, collection_name = 'hallmark', gene_set_name = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")

p <- via_sigs_plotting %>% ggplot(aes(`YM-155`, `BAX-channel-blocker`, color = "gray")) + 
  geom_point(size = 3) + 
  geom_point(data = via_sigs_plotting[via_sigs_plotting$gene_class == "cell cycle",], aes(`YM-155`, `BAX-channel-blocker`, color = "red"), size = 3) + 
  cdsr::theme_Publication() + 
  geom_label_repel(data = unique(rbind(
    via_sigs_plotting[via_sigs_plotting$gene_class == "cell cycle",] %>% arrange(-abs(`YM-155` * `BAX-channel-blocker`)) %>% head(10)
    # via_sigs_plotting[via_sigs_plotting$gene_class == "cell cycle",] %>% arrange(-abs(`BAX-channel-blocker`)) %>% head(5)
  )), aes(label = Gene, color = "red"), size = 8) + 
  # scale_color_manual(values=c("gray"="gray", "red"="red"), labels = c('c2','c1')) + 
  scale_colour_manual(name = 'Cell cycle-related', 
                      values =c('gray'='gray','red'='red'), labels = c('False', 'True')) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  l1000.analysis::theme_scatterplot() + 
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))
p
ggsave(file.path(plot_dir, "ym155_bax_comparison.png"), p, height = 10, width = 10)


p <- via_sigs_plotting %>% ggplot(aes(MCL1, NFKBIB, color = "gray")) + 
  geom_point(size = 3) + 
  geom_point(data = via_sigs_plotting[via_sigs_plotting$gene_class == "apoptosis",], aes(NFKBIB, MCL1, color = "red"), size = 3) + 
  cdsr::theme_Publication() + 
  geom_label_repel(data = unique(rbind(
    via_sigs_plotting[via_sigs_plotting$gene_class == "apoptosis",] %>% arrange(-abs(NFKBIB * MCL1)) %>% head(10)
    # via_sigs_plotting[via_sigs_plotting$gene_class == "cell cycle",] %>% arrange(-abs(`BAX-channel-blocker`)) %>% head(5)
  )), aes(label = Gene, color = "red"), size = 8) + 
  # scale_color_manual(values=c("gray"="gray", "red"="red"), labels = c('c2','c1')) + 
  scale_colour_manual(name = 'apoptosis-related', 
                      values =c('gray'='gray','red'='red'), labels = c('False', 'True')) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  l1000.analysis::theme_scatterplot() + 
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))
p
ggsave(file.path(plot_dir, "fgfr3_bcl2l1_comparison.png"), p, height = 10, width = 10)


##


### GSEA ON SINGLE CPDS ------

gsc_data <- l1000.analysis::load_gsc_data()
top_genes <- via_sigs_plotting$Gene[order(-via_sigs_plotting$`YM-155`)] %>% head(20)
gsa_out <- cdsr::run_GSAhyper(hit_genes = top_genes, all_genes = via_sigs_plotting$Gene, gsc = gsc_data[['hallmark']])
gsa_out %>% arrange(`p-value`) %>% head(10)


### LOOK AT VIABILITY DIFFERENCES BETWEEN THE CLUSTERS ------


trimmed_df <- merged_data[,c("auc_avg", "auc_prism", "auc_gdsc", "auc_ctrp", "tas", "cell_id", "pert_iname")] %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)
trimmed_df$cluster_assignment <- llply(trimmed_df$pert_iname, function(x) {
  if (x %in% weird_cluster_perts) {
    return("Cluster 2")
  } else if (x %in% main_cluster_perts) {
    return("Cluster 1")
  } else {
    return("")
  }
}) %>% as.character()
p <- trimmed_df %>% ggplot(aes(auc_avg, tas, color = cluster_assignment)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  scale_color_manual(values=c("gray", "red", "blue"))
ggExtra::ggMarginal(p, groupColour = T)







### GENE SET SCORES ------

gene_set_scores <- data.frame(
  proliferation = all_viability_sigs[proliferation_genes,] %>% colMeans(na.rm = T),
  apoptosis = all_viability_sigs[apoptosis_genes,] %>% colMeans(na.rm = T),
  cell_cycle = all_viability_sigs[cell_cycle_genes,] %>% colMeans(na.rm = T),
  pert_iname = colnames(all_viability_sigs)
)

gene_set_scores$cluster_assignment <- llply(gene_set_scores$pert_iname, function(x) {
  if (x %in% weird_cluster_perts) {
    return("Cluster 2")
  } else if (x %in% main_cluster_perts) {
    return("Cluster 1")
  } else {
    return("")
  }
}) %>% as.character()

p <- gene_set_scores %>% ggplot(aes(apoptosis, cell_cycle, color = cluster_assignment)) + 
  geom_point(size = 3) + 
  cdsr::theme_Publication() + 
  scale_color_manual(values=c("gray", "red", "blue")) + 
  geom_label_repel(data = unique(rbind(
    gene_set_scores %>% arrange(-abs(apoptosis)) %>% head(10),
    gene_set_scores %>% arrange(-abs(cell_cycle)) %>% head(10)#,
    # gene_set_scores[gene_set_scores$pert_iname %in% c(weird_cluster_perts, main_cluster_perts),]
  )), aes(label = pert_iname), size = 8) + 
  l1000.analysis::theme_scatterplot() + 
  labs(x = "Apoptosis NES", y = "Cell cycle NES", color = "Cluster") + 
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))
p
ggsave(file.path(plot_dir, "gene_set_avg_each_cluster.png"), p, height = 10, width = 10)








### LOOK AT TAS-VIABILITY RELATIONSHIPS ------

trimmed_data_collapsed <- merged_data[,c("cell_id", "pert_iname", "auc_avg", "tas")] %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

tas_via_cors <- ddply(.data = trimmed_data_collapsed, .variables = ~ pert_iname, .fun = function(x) {
  cor(x$auc_avg, x$tas)
}) %>% 
  set_colnames(c("pert_iname", "tas_via_cor"))

tas_via_cors$cluster_assignment <- llply(tas_via_cors$pert_iname, function(x) {
  if (x %in% weird_cluster_perts) {
    return("Cluster 2")
  } else if (x %in% main_cluster_perts) {
    return("Cluster 1")
  } else {
    return("")
  }
}) %>% as.character()


tas_via_cors_to_plot <- tas_via_cors[tas_via_cors$cluster_assignment != "",]
p <- tas_via_cors_to_plot %>% ggplot(aes(tas_via_cor, color = cluster_assignment)) + 
  geom_density(size = 3) + 
  cdsr::theme_Publication() + 
  scale_color_manual(name = 'Cell cycle-related', 
                      values = c('Cluster 1' = 'red', "Cluster 2" = "blue")) + 
  l1000.analysis::theme_scatterplot() +
  labs(x = "Corr. bt exp strength\nand sensivity") + 
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))
p
ggsave(file.path(plot_dir, "exp_viability_relationship_by_cluster.png"), p, height = 10, width = 10)



