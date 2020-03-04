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

# merged_data_collapsed <- merged_data %>% 
#   group_by(pert_iname, cell_id) %>% 
#   summarise_all(mean)
# merged_data_collapsed %<>% as.data.frame()




### COMPUTE VIABILITY SIGNATURES FOR EACH COMPOUND --------

cpd_viability_sigs <- l1000.analysis::compute_viability_signature_all_cpds(data_list = l1000_data, sens_column = "auc_avg")

cpd_viability_sigs_melted <- melt(cpd_viability_sigs)
viability_related_sigs <- cpd_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "t")
viability_related_sigs_cast <- reshape2::acast(viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()

pca_out <- prcomp(viability_related_sigs_cast)
pca_embed <- pca_out$x %>% as.data.frame() %>% rownames_to_column("gene")
pca_embed %>% 
  as.data.frame() %>% 
  ggplot(aes(PC1, PC2)) + 
  geom_point() + 
  geom_text_repel(data = unique(rbind(pca_embed %>% arrange(-abs(PC1)) %>% head(10),
                                      pca_embed %>% arrange(-abs(PC2)) %>% head(10))),
                  aes(label = gene))


l1000.analysis::make_via_sig_volcano(cpd_viability_sigs$penfluridol)
l1000.analysis::plot_gene_auc_relationship(cpd = "penfluridol", via_sig = cpd_viability_sigs$penfluridol, data_list = l1000_data, gene = "NOP56")


### HEATMAP VIA SIG CORRELATIONS -------

gene_info <- l1000.analysis::load_gene_info()
lm_genes <- gene_info$pr_gene_symbol[gene_info$pr_is_lm == 1]

via_sig_cors <- cor(viability_related_sigs_cast[lm_genes,])
diag(via_sig_cors) <- NA

pheatmap::pheatmap(via_sig_cors, 
                   filename = "~/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots/cpd_sig_heatmap.png", 
                   height = 5, 
                   width = 5, 
                   labels_row = F, labels_col = F, show_rownames = F, show_colnames = F)



corsello.drug.categories <- load.from.taiga(data.name='repurposing-compound-annotations-e8de', data.version=9, data.file='corsello_drug_categories')

annot_df <- merge(data.frame(row.names = rownames(via_sig_cors),
                 name = rownames(via_sig_cors)),
      corsello.drug.categories[,c("name", "drug_category", "is_cancer_merged")],
      all.x = T) %>% 
  column_to_rownames("name")


pheatmap::pheatmap(via_sig_cors, annotation_row = annot_df)



###

### LOOK AT RELATIONSHIP BETWEEN VIA SIG STRENGTH AND SELECTIVITY -------

via_sig_strengths <- l1000.analysis::compute_viability_signature_strengths(cpd_viability_sigs)
cpd_selectivities <- l1000.analysis::compute_sensitivity_selectivities(l1000_data)
cpd_avg_sens <- l1000.analysis::compute_avg_sensitivities(l1000_data)


sig_stats <- merge(via_sig_strengths, cpd_selectivities, by = "pert_iname")
sig_stats <- merge(sig_stats, cpd_avg_sens, by = "pert_iname")


sig_stats %>% ggplot(aes(selectivity, via_sig_strength, size = selectivity)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_label_repel(data = unique(rbind(sig_stats %>% arrange(-abs(via_sig_strength)) %>% head(10),
                                       sig_stats %>% arrange(-abs(selectivity)) %>% head(10),
                                       sig_stats[sig_stats$pert_iname == "nutlin-3",])),
                   aes(label = pert_iname))


tp53.mutation.data <- load.from.taiga(data.name='relevant-mutation-data-1424', data.file='tp53_mutation_data')
braf.mutation.data <- load.from.taiga(data.name='relevant-mutation-data-1424', data.file='braf_mutation_data')


cpd_of_interest <- "tanespimycin"
l1000.analysis::plot_pert_data_heatmap(data_list = l1000_data, via_sig_list = cpd_viability_sigs, cpd_name = cpd_of_interest) #, additional_metadata = braf.mutation.data)

l1000.analysis::make_via_sig_volcano(via_sig_output = cpd_viability_sigs[[cpd_of_interest]], cpd_name = cpd_of_interest)

##


### CALCULATE CORRELATIONS BETWEEN VIA-RELATED and -UNRELATED COMPONENTS --------


via_sig_cors <- ldply(cpd_viability_sigs, function(x) {
  merged_sigs <- merge(x$viability_related, x$viability_unrelated,
                       by = "Gene",
                       suffixes = c("_slope", "_intercept"))
  return(cor(merged_sigs$logFC_slope, merged_sigs$logFC_intercept))
}) %>% set_colnames(c("pert_iname", "via_sig_cors"))
via_sig_cors %>% arrange(abs(via_sig_cors)) %>% head(20)

##


### RUN GSEA ON A VIABILITY SIGNATURE ------

cpd_of_interest <- "YM-155"
viability_sig_of_interest <- cpd_viability_sigs[[cpd_of_interest]][["viability_related"]]

gsc_data <- l1000.analysis::load_gsc_data()
gsea_out <- cdsr::run_fGSEA(gsc = gsc_data[['hallmark']], gene_stat = viability_sig_of_interest$logFC %>% set_names(viability_sig_of_interest$Gene), perm_type = "gene")
gsea_out %>% arrange(-abs(NES)) %>% head(10)




### PLOT EXAMPLE EXPRESSION-AUC RELATIONSHIP -------


cpd_of_interest <- "PD-98059"
p <- l1000.analysis::plot_gene_auc_relationship(cpd = cpd_of_interest, 
                                           via_sig = cpd_viability_sigs[[cpd_of_interest]], 
                                           data_list = l1000_data, 
                                           gene = "RRM2", 
                                           sens_column = "auc_avg", is_plot_abline = F, is_plot_intercept_dot = T) + 
  labs(y = "Gene DE")
p

ggsave(file.path(plot_dir, "example_gene_auc_relationship.png"), p, height = 10, width = 11)










### PLOT VOLCANO PLOTS OF VIABILITY-RELATED AND -UNRELATED COMPONENTS --------

cpd_to_plot <- "tanespimycin"

p <- l1000.analysis::make_via_sig_volcano(via_sig_output = cpd_viability_sigs[[cpd_to_plot]], cpd_name = cpd_to_plot)
p
ggsave(file.path(plot_dir, "nutlin_volcano_sigs.png"), p, height = 10, width = 20)

l1000.analysis::plot_via_sig_comparison(via_sig_output = cpd_viability_sigs[[cpd_to_plot]])

gsc_data <- l1000.analysis::load_gsc_data()
e2f_genes <- l1000.analysis::get_gene_set_genes(gsc = gsc_data, collection_name = "hallmark", gene_set_name = "HALLMARK_E2F_TARGETS")


### PLOT VIABILITY-RELATED AND -UNRELATED AGAINST EACH OTHER --------

l1000.analysis::plot_via_sig_comparison(via_sig_output = cpd_viability_sigs$`PD-98059`)



cpd_to_plot <- "tanespimycin"
gene_class_df <- ldply(cpd_viability_sigs$tanespimycin$viability_unrelated$Gene, function(x) {
  if (x %in% c("HSPA8", "HSPA6", "HSPD1", "DNAJB1", "HSPA1A")) {
    return("MOA-related")
  } else if (x %in% c("PTK2", "MAL", "PSMD2")) {
    return("Death-related")
  } else {
    return("")
  }
}) %>%
  set_colnames("gene_class") %>%
  mutate(Gene = cpd_viability_sigs$tanespimycin$viability_unrelated$Gene)


p <- l1000.analysis::plot_via_sig_comparison(via_sig_output = cpd_viability_sigs[[cpd_to_plot]], gene_class_df = gene_class_df, cpd_to_plot)
ggsave(file.path(plot_dir, "via_sig_comparison_tanespimycin.png"), p, height = 10, width = 10)


cpd_to_plot <- "PD-98059"
gene_class_df <- ldply(cpd_viability_sigs$tanespimycin$viability_unrelated$Gene, function(x) {
  if (x %in% c("DUSP4", "DUSP6", "TIPARP", "IER3", "FOSL1")) {
    return("MOA-related")
  } else if (x %in% c("PUF60", "CLTC")) {
    return("Death-related")
  } else {
    return("")
  }
}) %>% 
  set_colnames("gene_class") %>% 
  mutate(Gene = cpd_viability_sigs$tanespimycin$viability_unrelated$Gene)


p <- l1000.analysis::plot_via_sig_comparison(via_sig_output = cpd_viability_sigs[[cpd_to_plot]], gene_class_df = gene_class_df, cpd_to_plot)
ggsave(file.path(plot_dir, "via_sig_comparison_other_drug.png"), p, height = 10, width = 10)


### PLOT GENE-AUC RELATIONSHIPS -------

cpd_to_plot <- "tanespimycin"

p1 <- l1000.analysis::plot_gene_auc_relationship(cpd = cpd_to_plot, via_sig = cpd_viability_sigs[[cpd_to_plot]], data_list = l1000_data, gene = "PTK2", sens_column = "auc_avg", 
                                           is_plot_intercept_dot = F, is_plot_abline = F)

p2 <- l1000.analysis::plot_gene_auc_relationship(cpd = cpd_to_plot, via_sig = cpd_viability_sigs[[cpd_to_plot]], data_list = l1000_data, gene = "HSPA6", sens_column = "auc_avg", 
                                           is_plot_intercept_dot = F, is_plot_abline = F)


p_grid <- cowplot::plot_grid(p1, p2)

ggsave(file.path(plot_dir, "gene_auc_relationships_drug.png"), p_grid, height = 10, width = 20)










### LOOK AT EACH CPD'S RELATIONSHIP TO GLOBAL VIA SIG --------


global_viability_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_avg")
global_via_related_sig <- global_viability_sig$viability_related[,c("Gene", "logFC")] %>% 
  column_to_rownames("Gene")

shared_genes <- intersect(rownames(global_via_related_sig), rownames(viability_related_sigs_cast))


global_cpd_cors <- cor(viability_related_sigs_cast[shared_genes,], global_via_related_sig[shared_genes,]) %>% 
  as.data.frame() %>% 
  set_colnames(c("global_cor"))


global_cpd_cors_annot <- merge(global_cpd_cors, repurp_data[,c("Name", "MOA")], by.x = "row.names", by.y = "Name", all.x = T)
global_cpd_cors_annot <- merge(global_cpd_cors_annot, corsello.drug.categories, by.x = "Row.names", by.y = "name", all.x = T)

global_cpd_cors_annot$MOA <- substr(global_cpd_cors_annot$MOA, 1, 20)


global_cpd_cors_annot %>% ggplot(aes(reorder(MOA, -global_cor), global_cor)) + 
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_point(aes(fill = drug_category), size = 5, shape = 21, position = position_jitterdodge()) + 
  labs(x = "", y = "Model Pearson correlation")






### DEFINE "CHEMO" VIA SIG AND COMPARE TO OTHERS ------

cpd_viability_sigs_melted <- melt(cpd_viability_sigs)
viability_related_sigs <- cpd_viability_sigs_melted %>% 
  filter(L2 == "viability_related" & variable == "t")
viability_related_sigs_cast <- reshape2::acast(viability_related_sigs, formula = Gene ~ L1) %>% 
  as.data.frame()



corsello.drug.categories <- load.from.taiga(data.name='repurposing-compound-annotations-e8de', data.version=9, data.file='corsello_drug_categories')
chemo_drugs <- corsello.drug.categories$name[corsello.drug.categories$drug_category == "chemo"]

chemo_idx <- which((l1000_data$metadata$pert_iname %in% chemo_drugs) == T)
l1000_data_chemo <- list(
  expression = l1000_data$expression[chemo_idx,],
  sensitivity = l1000_data$sensitivity[chemo_idx,],
  metadata = l1000_data$metadata[chemo_idx,]
)


global_viability_sig_chemo <- l1000.analysis::compute_viability_signature(data_list = l1000_data_chemo, global = T, sens_column = "auc_avg")

via_related_sig_chemo <- global_viability_sig_chemo$viability_related %>% 
  set_rownames(global_viability_sig_chemo$viability_related$Gene)


chemo_via_sig_cors <- cor(viability_related_sigs_cast[rownames(via_related_sig_chemo),], via_related_sig_chemo$t) %>% 
  as.data.frame() %>% 
  set_colnames("chemo_cor") %>% 
  rownames_to_column("pert_iname") %>% 
  left_join(corsello.drug.categories[,c("name", "moa", "target", "drug_category")], by = c("pert_iname" = "name"))

pert_metrics <- load.from.taiga(data.name='l1000-metadata-8489', data.version=1, data.file='GSE92742_Broad_LINCS_pert_metrics')

chemo_via_sig_cors %<>% left_join(pert_metrics, by = c("pert_iname"))


chemo_via_sig_cors %>% arrange(chemo_cor)

chemo_via_sig_cors$moa <- substr(chemo_via_sig_cors$moa, 1, 20)
chemo_via_sig_cors %>% ggplot(aes(reorder(moa, -chemo_cor), chemo_cor)) + 
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_point(aes(fill = drug_category), size = 5, shape = 21, position = position_jitterdodge()) + 
  labs(x = "", y = "Model Pearson correlation")


chemo_via_sig_cors %>% ggplot(aes(tas_q75, chemo_cor)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  cdsr::theme_Publication()



via_sig_cors <- cor(viability_related_sigs_cast[lm_genes,])
diag(via_sig_cors) <- NA

chemo_via_sig_cors <- chemo_via_sig_cors[!duplicated(chemo_via_sig_cors$pert_iname),]
chemo_via_sig_cors %<>% set_rownames(chemo_via_sig_cors$pert_iname)
annot_df <- chemo_via_sig_cors[rownames(via_sig_cors), c("chemo_cor", "tas_q75"), drop = F]

pheatmap::pheatmap(via_sig_cors, annotation_row = annot_df, fontsize = 5)

pheatmap::pheatmap(via_sig_cors, annotation_row = annot_df, fontsize = 5, height = 10, width = 10, filename = "~/Desktop/a.png")



tsne_out <- Rtsne::Rtsne(viability_related_sigs_cast[lm_genes,] %>% t())

stopifnot(all.equal(colnames(viability_related_sigs_cast), rownames(chemo_via_sig_cors)))

tsne_coords <- tsne_out$Y %>% 
  as.data.frame() %>% 
  set_colnames(c("t1", "t2")) %>% 
  mutate(pert_iname = colnames(viability_related_sigs_cast)) %>% 
  left_join(chemo_via_sig_cors, by = "pert_iname")
  
tsne_coords %>% ggplot(aes(t1, t2, color = tas_q75)) + 
  geom_point(size = 3) + 
  cdsr::theme_Publication()









### LOOK AT DIFFERENCE BETWEEN (CUSTOM) CLUSTER SIGS ------

cluster1_cpds <- c("teniposide", "YM-155", "PF-04217903", "simvastatin")
cluster2_cpds <- c("flutamide", "rucaparib", "BAX-channel-blocker", "sonidegib")

sig_df <- data.frame(
  cluster1_sig = rowMeans(viability_related_sigs_cast[,cluster1_cpds]),
  cluster2_sig = rowMeans(viability_related_sigs_cast[,cluster2_cpds])
) %>% 
  rownames_to_column("Gene")

sig_df %>% ggplot(aes(cluster1_sig, cluster2_sig)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  geom_label_repel(data = unique(rbind(sig_df %>% arrange(-abs(cluster1_sig)) %>% head(20),
                                       sig_df %>% arrange(-abs(cluster2_sig)) %>% head(20))),
                   aes(label = Gene))


# gsc_data <- l1000.analysis::load_gsc_data()
# gsea_out <- cdsr::run_fGSEA(gene_stat = sig_df$cluster1_sig %>% set_names(sig_df$Gene), gsc = gsc_data[['hallmark']], perm_type = "gene")
# gsea_out %>% arrange(-NES) %>% head(20)



# 
# apoptosis_genes <- read.table("~/Documents/cancer_data_science/one_offs/uri_analysis/data/apoptosis_genes.txt")
# cc_genes <- read.table("~/Documents/cancer_data_science/one_offs/uri_analysis/data/cell_cycle_genes.txt")
# proliferation_genes <- read.table("~/Documents/cancer_data_science/one_offs/uri_analysis/data/proliferation_genes.txt")
# sig_df$cluster1_sig[sig_df$Gene %in% apoptosis_genes$V2] %>% mean()
# sig_df$cluster2_sig[sig_df$Gene %in% apoptosis_genes$V2] %>% mean()
# sig_df$cluster2_sig[sig_df$Gene %in% proliferation_genes$V2] %>% mean()
# 
# handcrafted_gsea <- data.frame(
#   apoptosis = viability_related_sigs_cast[rownames(viability_related_sigs_cast) %in% apoptosis_genes$V2,] %>% 
#     colMeans(),
#   cell_cycle = viability_related_sigs_cast[rownames(viability_related_sigs_cast) %in% cc_genes$V2,] %>% 
#     colMeans(),
#   proliferation = viability_related_sigs_cast[rownames(viability_related_sigs_cast) %in% proliferation_genes$V2,] %>% 
#     colMeans(),
#   pert_iname = colnames(viability_related_sigs_cast)
# )
# handcrafted_gsea$reverse_sig = handcrafted_gsea$pert_iname %in% cluster2_cpds
# handcrafted_gsea %<>% left_join(corsello.drug.categories, by = c("pert_iname" = "name"))
# 
# p <- handcrafted_gsea %>% ggplot(aes(apoptosis, proliferation, color = reverse_sig)) + 
#   geom_point(size = 3) + 
#   cdsr::theme_Publication() + 
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) + 
#   geom_label_repel(data = unique(rbind(handcrafted_gsea %>% arrange(-abs(apoptosis)) %>% head(10),
#                                        handcrafted_gsea %>% arrange(-abs(proliferation)) %>% head(10),
#                                        handcrafted_gsea %>% filter(reverse_sig == T))),
#                                        # handcrafted_gsea %>% filter(drug_category == "chemo"))), 
#                    aes(label = pert_iname),
#                    size = 6) +
#   l1000.analysis::theme_scatterplot() + 
#   theme(legend.text=element_text(size=20), legend.title=element_text(size=20)) + 
#   labs(x = "Apoptosis", y = "Cell cycle", color = '"Proliferative" cluster')
# p
# ggsave("~/Desktop/manual_gsea.png", p, height = 10, width = 10)
# 
# 
# l1000_data$sensitivity[l1000_data$metadata$pert_iname %in% cluster2_cpds,] %>% 
#   mutate(pert_iname = l1000_data$metadata$pert_iname[l1000_data$metadata$pert_iname %in% cluster2_cpds]) %>% 
#   arrange(pert_iname)
# 
# 
# primary.merged.MEDIAN.LFCVC.CB <- load.from.taiga(data.name='primary-screen-acc3', data.version=9, data.file='primary_merged_MEDIAN_LFCVC_CB') %>% t() %>% 
#   as.data.frame()
# primary.merged.collapse.col.meta <- load.from.taiga(data.name='primary-screen-acc3', data.version=9, data.file='primary_merged_collapse_col_meta') %>% as.data.frame()
# primary.merged.row.meta <- load.from.taiga(data.name='primary-screen-acc3', data.version=9, data.file='primary_merged_row_meta')
# 
# stopifnot(all.equal(rownames(primary.merged.MEDIAN.LFCVC.CB), primary.merged.collapse.col.meta$profile_id))
# 
# primary.merged.MEDIAN.LFCVC.CB %<>% set_rownames(primary.merged.collapse.col.meta$profile_id)
# 
# colnames(primary.merged.MEDIAN.LFCVC.CB) <- llply(colnames(primary.merged.MEDIAN.LFCVC.CB), function(x) {
#   curr_c <- primary.merged.row.meta$ccle_name[primary.merged.row.meta$feature_id == x]
#   return(strsplit(curr_c, "_")[[1]][1])
# }) %>% as.character()
# 
# 
# drug_to_get <- "sonidegib"
# drug_data <- primary.merged.MEDIAN.LFCVC.CB[primary.merged.collapse.col.meta$repurposing_name == drug_to_get,]
# drug_data <- drug_data[rowSums(drug_data %>% is.na()) < dim(drug_data)[2],]
# drug_data_melt <- drug_data %>% reshape2::melt() %>% 
#   set_colnames(c("cell_id", "dose_level_viability"))
# drug_metadata <- primary.merged.collapse.col.meta[primary.merged.collapse.col.meta$repurposing_name == drug_to_get,]
# 
# drug_l1k_data <- merged_data[merged_data$pert_iname == drug_to_get,] %>% 
#   group_by(cell_id) %>% 
#   summarise_all(mean)
# 
# drug_l1k_data_with_dose_via <- merge(drug_l1k_data, drug_data_melt, by = "cell_id")
# drug_l1k_data_with_dose_via %>% ggplot(aes(dose_level_viability, TUBA4A)) + 
#   geom_point() + 
#   cdsr::theme_Publication()
# 
# drug_l1k_data_with_dose_via %>% ggplot(aes(dose_level_viability, auc_ctrp)) + 
#   geom_point() + 
#   cdsr::theme_Publication()
# 
# 
# 
# 
# 
# gene_info <- l1000.analysis::load_gene_info()
# lm_genes <- gene_info$pr_gene_symbol[gene_info$pr_is_lm == 1]
# viability_related_sigs_cast_plotting <- viability_related_sigs_cast %>% 
#   rownames_to_column("Gene")
# viability_related_sigs_cast_plotting %>% ggplot(aes(vemurafenib, teniposide, color = Gene %in% lm_genes)) + 
#   geom_point() + 
#   cdsr::theme_Publication() + 
#   geom_label_repel(data = unique(rbind(viability_related_sigs_cast_plotting %>% arrange(-abs(vemurafenib)) %>% head(20),
#                                        viability_related_sigs_cast_plotting %>% arrange(-abs(teniposide)) %>% head(20))),
#                    aes(label = Gene)) + 
#   geom_smooth(method = "lm")
# 
# 
# l1000.analysis::plot_gene_auc_relationship(cpd = "rucaparib", via_sig = cpd_viability_sigs$rucaparib, data_list = l1000_data, gene = "TUBA4A", sens_column = "auc_avg")
# 
# merged_data[merged_data$pert_iname == "rucaparib", c("auc_avg", "auc_prism", "auc_gdsc", "auc_ctrp", "cell_id")] %>% arrange(auc_avg)
# 
# masterfile.2019.06.14 <- load.from.taiga(data.name='master-cell-line-export-0306', data.version=339, data.file='masterfile_2019-06-14')
# masterfile.2019.06.14$cell_id <- llply(masterfile.2019.06.14$CCLE_name, function(x) {
#   strsplit(x, "_")[[1]][1]
# }) %>% as.character()
# masterfile.2019.06.14 <- masterfile.2019.06.14[!duplicated(masterfile.2019.06.14$cell_id),]
# 
# data_of_interest <- merged_data[merged_data$pert_iname == "flutamide", c("auc_avg", "auc_prism", "auc_gdsc", "auc_ctrp", "cell_id", "TUBA4A")] %>% 
#   group_by(cell_id) %>% 
#   summarise_all(mean) %>% 
#   left_join(masterfile.2019.06.14[,c("cell_id", "Primary Disease")], by = "cell_id")
#   
# data_of_interest %>% ggplot(aes(auc_avg, TUBA4A, color = `Primary Disease`)) + 
#   geom_point(size = 3) + 
#   cdsr::theme_Publication() + 
#   geom_label_repel(data = data_of_interest, aes(label = `Primary Disease`)) + 
#   geom_label_repel(data = data_of_interest, aes(label = cell_id))
# 
# 
# data_of_interest_full <- merged_data[merged_data$pert_iname == "teniposide",] %>% 
#   group_by(cell_id) %>% 
#   summarise_all(mean) %>% 
#   left_join(masterfile.2019.06.14[,c("cell_id", "Primary Disease")], by = "cell_id")
# cor(data_of_interest_full[data_of_interest_full$cell_id == "THP1", rownames(viability_related_sigs_cast)] %>% t(), viability_related_sigs_cast$teniposide)
# plot(data_of_interest_full[data_of_interest_full$cell_id == "THP1", rownames(viability_related_sigs_cast)] %>% t(), viability_related_sigs_cast$teniposide)
# 
# 
# cor(merged_data[(merged_data$pert_iname == "rucaparib") & (merged_data$cell_id %in% c("SNGM", "SNUC5")), rownames(viability_related_sigs_cast)] %>% colMeans(),
#     merged_data[(merged_data$pert_iname == "flutamide") & (merged_data$cell_id == "CL34"), rownames(viability_related_sigs_cast)] %>% t())
# 
# 
# 
# 
# gsc_data <- l1000.analysis::load_gsc_data()
# gsea_out <- cdsr::run_fGSEA(gene_stat = viability_related_sigs_cast$flutamide %>% set_names(viability_related_sigs_cast %>% rownames()), gsc = gsc_data[['hallmark']], perm_type = "gene")
# gsea_out %>% arrange(-abs(NES)) %>% head(10)
# 
# 
# weird_cpds <- c("flutamide", "rucaparib", "BAX-channel-blocker", "sonidegib")
# top_genes <- rownames(viability_related_sigs_cast)[order(-abs(rowMeans(viability_related_sigs_cast[,weird_cpds])))[1:30]]
# pheatmap::pheatmap(scale(viability_related_sigs_cast)[top_genes, weird_cpds])



