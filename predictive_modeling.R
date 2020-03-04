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
library(glmnet)

plot_dir <- "~/Documents/cancerdatascience/l1000_via_sig_ms/figures/plots"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data <- cbind(l1000_data$expression, l1000_data$sensitivity, l1000_data$metadata)

gene_info <- l1000.analysis::load_gene_info()
lm_genes <- gene_info$pr_gene_symbol[gene_info$pr_is_lm == 1]

merged_data_collapsed <- merged_data[,c(lm_genes, "auc_avg", "pert_iname", "cell_id", "tas")] %>% 
  group_by(pert_iname, cell_id) %>% 
  summarise_all(mean)

merged_data_collapsed %<>% as.data.frame()


### CHECK NUMBER OF CELL LINES FOR EACH COMPOUND -----

pert_counts <- table(merged_data_collapsed$pert_iname) %>% 
  as.data.frame() %>% 
  set_colnames(c("pert_iname", "count"))

min(pert_counts$count)






### RUN ELASTIC NET ON GLOBAL MODEL --------




n_fold_cv <- 10
folds <- cvTools::cvFolds(n = length(unique(merged_data_collapsed$pert_iname)), K = n_fold_cv)
unique_perts <- unique(merged_data_collapsed$pert_iname)
lambda_val <- 0.003
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

for (fold_ii in seq(n_fold_cv)) {
  
  ## split data
  curr_fold_test_idx <- which(folds$which == fold_ii)
  
  curr_train_perts <- unique_perts[-curr_fold_test_idx]
  curr_test_perts <- unique_perts[curr_fold_test_idx]
  
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, lm_genes] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, "auc_avg"] %>% as.matrix()
  
  curr_test_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, lm_genes] %>% as.matrix()
  curr_test_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
}

lm_pred_df <- data.frame(
  lm_type = c(
    rep("expression_kfold", length(truth_pred_cor_list)),
    rep("tas_kfold", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  )
)

lm_pred_df %>% ggplot(aes(lm_type, results)) + 
  geom_boxplot() +
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_point(size = 5, shape = 21, position = position_jitter()) + 
  labs(x = "", y = "Model R")


lm_pred_df[lm_pred_df$lm_type == "expression_kfold", "results"] %>% mean()







### RUN ELASTIC NET PER-COMPOUND -------


per_cpd_variability <- merged_data_collapsed[,c("pert_iname", "auc_avg")] %>% 
  group_by(pert_iname) %>% 
  summarise_all(sd) %>% 
  as.data.frame() %>% 
  set_colnames(c("cpd", "sd"))

hist(per_cpd_variability$sd, 20)

sd_cutoff <- 0.0
perts_with_sufficient_sd <- per_cpd_variability$cpd[per_cpd_variability$sd >= sd_cutoff]

data_for_lm <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% perts_with_sufficient_sd,]


unique_perts <- unique(data_for_lm$pert_iname)
lambda_val <- 0.006
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

pb = txtProgressBar(min = 0, max = length(unique_perts), initial = 0) 

for (cpd_ii in seq(length(unique_perts))) {
  
  
  ## split data
  curr_test_cpd <- unique_perts[cpd_ii]
  curr_train_cpds <- unique_perts[-which(unique_perts == curr_test_cpd)]
  
  # use full data for training
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, lm_genes] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, "auc_avg"] %>% as.matrix()
  
  # use only selective compounds for testing
  curr_test_data <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, lm_genes] %>% as.matrix()
  curr_test_labels <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = data_for_lm[data_for_lm$pert_iname %in% curr_train_cpds, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
  setTxtProgressBar(pb, cpd_ii)
  
  
}

lm_pred_per_cpd_df <- data.frame(
  lm_type = c(
    rep("Expression", length(truth_pred_cor_list)),
    rep("TAS", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  ),
  cpd = c(rep(unique_perts, 2))
)

lm_pred_per_cpd_df %>% ggplot(aes(lm_type, results)) + 
  geom_boxplot() +
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_point(size = 5, shape = 21, position = position_jitter()) + 
  labs(x = "", y = "Model R")


lm_pred_per_cpd_df_cast <- reshape2::acast(data = lm_pred_per_cpd_df, formula = cpd ~ lm_type, value.var = "results", fun.aggregate = mean) %>% 
  as.data.frame() %>% 
  rownames_to_column("cpd")


lm_pred_per_cpd_df_cast %>% ggplot(aes(Expression, TAS)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_abline()



### PLOT BOTH RESULTS TOGETHER -------


combined_results <- rbind(lm_pred_df, lm_pred_per_cpd_df[,1:2])


p <- combined_results %>% ggplot(aes(lm_type, results)) +
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  geom_jitter(shape=16, size = 4, position=position_jitter(0.2)) + 
  scale_x_discrete(labels=c("expression_kfold" = "Expression",
                            "tas_kfold" = "TAS",
                            "Expression" = "Expression",
                            "TAS" = "TAS")) + 
  labs(x = "", y = "Model Pearson corr.") + 
  l1000.analysis::theme_scatterplot() + 
  geom_vline(xintercept=c(2.5))
p
ggsave(file.path(plot_dir, "predictive_model_boxplot2.png"), p, height = 10, width = 13)




### LOOK AT BREAKDOWN BY MOA -----


repurp_data <- l1000.analysis::load_repurposing_cpd_annotations()
repurp_data <- repurp_data[!duplicated(repurp_data$name),]
lm_pred_per_cpd_df_annot <- merge(repurp_data[,c("name", "moa")], lm_pred_per_cpd_df, by.x = "name", by.y = "cpd", all.y = T)
lm_pred_per_cpd_df_annot <- lm_pred_per_cpd_df_annot[lm_pred_per_cpd_df_annot$lm_type == "Expression",]

moa_counts <- table(lm_pred_per_cpd_df_annot$moa)
top_moas <- names(moa_counts)[moa_counts > 0]
lm_pred_per_cpd_df_annot_common <- lm_pred_per_cpd_df_annot[lm_pred_per_cpd_df_annot$moa %in% top_moas,]

corsello.drug.categories <- load.from.taiga(data.name='repurposing-compound-annotations-e8de', data.version=9, data.file='corsello_drug_categories')

lm_pred_per_cpd_df_annot_common <- merge(lm_pred_per_cpd_df_annot_common, corsello.drug.categories[,c("drug_category", "name", "is_cancer_merged", "disease.area")], by.x = "name", by.y = "name", all.x = T)

lm_pred_per_cpd_df_annot_common$moa <- substr(lm_pred_per_cpd_df_annot_common$moa, 1, 20)

lm_pred_per_cpd_df_annot_common %>% ggplot(aes(reorder(moa, -results), results)) + 
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_point(aes(fill = drug_category), size = 5, shape = 21, position = position_jitterdodge()) + 
  labs(x = "", y = "Model Pearson correlation")

lm_pred_per_cpd_df_annot_common %>% ggplot(aes(drug_category, results)) + 
  geom_boxplot() + 
  cdsr::theme_Publication() + 
  geom_point(aes(fill = drug_category), size = 5, shape = 21, position = position_jitterdodge())
  


##



### TRY PREDICTING WITH SUBSETS OF DRUG CATEGORIES ------

# no chemo


per_cpd_variability <- merged_data_collapsed[,c("pert_iname", "auc_avg")] %>% 
  group_by(pert_iname) %>% 
  summarise_all(sd) %>% 
  as.data.frame() %>% 
  set_colnames(c("cpd", "sd"))

hist(per_cpd_variability$sd, 20)

sd_cutoff <- 0.10
perts_with_sufficient_sd <- per_cpd_variability$cpd[per_cpd_variability$sd >= sd_cutoff]

data_for_lm <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% perts_with_sufficient_sd,]

chemo_drugs <- corsello.drug.categories$name[corsello.drug.categories$drug_category == "chemo"]
data_for_lm <- data_for_lm[!(data_for_lm$pert_iname %in% chemo_drugs),]


unique_perts <- unique(data_for_lm$pert_iname)
lambda_val <- 0.006
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

pb = txtProgressBar(min = 0, max = length(unique_perts), initial = 0) 

for (cpd_ii in seq(length(unique_perts))) {
  
  
  ## split data
  curr_test_cpd <- unique_perts[cpd_ii]
  curr_train_cpds <- unique_perts[-which(unique_perts == curr_test_cpd)]
  
  # use full data for training
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, lm_genes] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, "auc_avg"] %>% as.matrix()
  
  # use only selective compounds for testing
  curr_test_data <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, lm_genes] %>% as.matrix()
  curr_test_labels <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = data_for_lm[data_for_lm$pert_iname %in% curr_train_cpds, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
  setTxtProgressBar(pb, cpd_ii)
  
  
}

lm_pred_per_cpd_no_chemo_df <- data.frame(
  lm_type = c(
    rep("Expression", length(truth_pred_cor_list)),
    rep("TAS", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  ),
  cpd = c(rep(unique_perts, 2))
)

lm_pred_per_cpd_no_chemo_df %>% ggplot(aes(lm_type, results)) + 
  geom_boxplot() +
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  geom_point(size = 5, shape = 21, position = position_jitter()) + 
  labs(x = "", y = "Model R")

both_results <- merge(lm_pred_per_cpd_df[lm_pred_per_cpd_df$lm_type == "Expression",], 
      lm_pred_per_cpd_no_chemo_df[lm_pred_per_cpd_no_chemo_df$lm_type == "Expression",], 
      by = "cpd",
      suffixes = c("_all", "_nochemo"))

both_results %>% ggplot(aes(results_all, results_nochemo)) + 
  geom_point() + 
  geom_abline() + 
  geom_text_repel(data = both_results %>% arrange(-abs(results_all - results_nochemo)) %>% head(10),
                  aes(label = cpd)) 


##


### PREDICT USING JUST APOPTOSIS AND CELL CYCLE GENES -------


per_cpd_variability <- merged_data_collapsed[,c("pert_iname", "auc_avg")] %>% 
  group_by(pert_iname) %>% 
  summarise_all(sd) %>% 
  as.data.frame() %>% 
  set_colnames(c("cpd", "sd"))

hist(per_cpd_variability$sd, 20)

sd_cutoff <- 0.10
perts_with_sufficient_sd <- per_cpd_variability$cpd[per_cpd_variability$sd >= sd_cutoff]

data_for_lm <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% perts_with_sufficient_sd,]

gsc <- l1000.analysis::load_gsc_data()
cell_cycle_genes <- c(l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_G2M_CHECKPOINT"),
                      l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_E2F_TARGETS"))
apoptosis_genes <- l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_APOPTOSIS")

genes_of_interest <- intersect(colnames(merged_data_collapsed), c(apoptosis_genes, cell_cycle_genes))




unique_perts <- unique(data_for_lm$pert_iname)
lambda_val <- 0.006
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

pb = txtProgressBar(min = 0, max = length(unique_perts), initial = 0) 

for (cpd_ii in seq(length(unique_perts))) {
  
  
  ## split data
  curr_test_cpd <- unique_perts[cpd_ii]
  curr_train_cpds <- unique_perts[-which(unique_perts == curr_test_cpd)]
  
  # use full data for training
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, genes_of_interest] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, "auc_avg"] %>% as.matrix()
  
  # use only selective compounds for testing
  curr_test_data <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, genes_of_interest] %>% as.matrix()
  curr_test_labels <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = data_for_lm[data_for_lm$pert_iname %in% curr_train_cpds, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
  setTxtProgressBar(pb, cpd_ii)
  
  
}

lm_pred_per_cpd_specific_genes_df <- data.frame(
  lm_type = c(
    rep("Expression", length(truth_pred_cor_list)),
    rep("TAS", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  ),
  cpd = c(rep(unique_perts, 2))
)


both_results <- merge(lm_pred_per_cpd_df %>% filter(lm_type == "Expression"), 
                      lm_pred_per_cpd_specific_genes_df %>% filter(lm_type == "Expression"), 
                      by = "cpd", 
                      suffixes = c("_lmgenes", "_specificgenes"))

both_results_annot <- merge(both_results, corsello.drug.categories, by.x = "cpd", by.y = "name", all.x = T)

both_results_annot %>% ggplot(aes(results_lmgenes, results_specificgenes, color = drug_category %>% as.factor())) + 
  geom_point(size = 3) + 
  cdsr::theme_Publication() + 
  geom_abline() + 
  l1000.analysis::theme_scatterplot() + 
  labs(x = "LM genes", y = "Curated gene sets") + 
  geom_text_repel(
    # data = unique(rbind(both_results %>% arrange(-abs(results_lmgenes)) %>% head(10),
    #                                   both_results %>% arrange(-abs(results_specificgenes)) %>% head(10))),
    data = both_results_annot %>% arrange(-abs(results_lmgenes - results_specificgenes)) %>% head(20),
                  aes(label = cpd))


cpd_viability_related_sigs_cast %>% 
  ggplot(aes(mesoridazine, flutamide, color = rownames(cpd_viability_related_sigs_cast) %in% genes_of_interest)) + 
  geom_point() + 
  cdsr::theme_Publication()


##




### PREDICT USING JUST TOP GLOBAL VIA SIG GENES -------


per_cpd_variability <- merged_data_collapsed[,c("pert_iname", "auc_avg")] %>% 
  group_by(pert_iname) %>% 
  summarise_all(sd) %>% 
  as.data.frame() %>% 
  set_colnames(c("cpd", "sd"))

hist(per_cpd_variability$sd, 20)

sd_cutoff <- 0.10
perts_with_sufficient_sd <- per_cpd_variability$cpd[per_cpd_variability$sd >= sd_cutoff]

data_for_lm <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% perts_with_sufficient_sd,]


global_viability_sig <- l1000.analysis::compute_viability_signature(data_list = l1000_data, global = T, sens_column = "auc_avg")
num_top_genes <- 20
top_sig_genes <- global_viability_sig$viability_related$Gene[order(-abs(global_viability_sig$viability_related$logFC))][1:num_top_genes]

genes_of_interest <- intersect(lm_genes, top_sig_genes)#, apoptosis_genes))




unique_perts <- unique(data_for_lm$pert_iname)
lambda_val <- 0.006
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

pb = txtProgressBar(min = 0, max = length(unique_perts), initial = 0) 

for (cpd_ii in seq(length(unique_perts))) {
  
  
  ## split data
  curr_test_cpd <- unique_perts[cpd_ii]
  curr_train_cpds <- unique_perts[-which(unique_perts == curr_test_cpd)]
  
  # use full data for training
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, genes_of_interest] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_cpds, "auc_avg"] %>% as.matrix()
  
  # use only selective compounds for testing
  curr_test_data <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, genes_of_interest] %>% as.matrix()
  curr_test_labels <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = data_for_lm[data_for_lm$pert_iname %in% curr_train_cpds, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
  setTxtProgressBar(pb, cpd_ii)
  
  
}

lm_pred_per_cpd_specific_genes_df <- data.frame(
  lm_type = c(
    rep("Expression", length(truth_pred_cor_list)),
    rep("TAS", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  ),
  cpd = c(rep(unique_perts, 2))
)


both_results <- merge(lm_pred_per_cpd_df %>% filter(lm_type == "Expression"), 
                      lm_pred_per_cpd_specific_genes_df %>% filter(lm_type == "Expression"), 
                      by = "cpd", 
                      suffixes = c("_lmgenes", "_specificgenes"))

both_results %>% ggplot(aes(results_lmgenes, results_specificgenes)) + 
  geom_point() + 
  cdsr::theme_Publication() + 
  geom_abline() + 
  l1000.analysis::theme_scatterplot() + 
  labs(x = "LM genes", y = "Curated gene sets") + 
  geom_text_repel(
    # data = unique(rbind(both_results %>% arrange(-abs(results_lmgenes)) %>% head(10),
    #                                   both_results %>% arrange(-abs(results_specificgenes)) %>% head(10))),
    data = both_results %>% arrange(-abs(results_lmgenes - results_specificgenes)) %>% head(20),
    aes(label = cpd))

results_melted <- melt(both_results)
results_melted %>% ggplot(aes(variable, value)) +
  geom_boxplot() + 
  cdsr::theme_Publication()


##



### GLOBAL MODEL ON TOP GLOBAL GENES -----

num_top_genes <- 200
lm_via_sig <- global_viability_sig$viability_related[global_viability_sig$viability_related$Gene %in% lm_genes,]
top_sig_genes <- lm_via_sig$Gene[order(-abs(lm_via_sig$logFC))][1:num_top_genes]

# genes_of_interest <- intersect(lm_genes, top_sig_genes)#, apoptosis_genes))

gsc <- l1000.analysis::load_gsc_data()
cell_cycle_genes <- c(l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_G2M_CHECKPOINT"),
                      l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_E2F_TARGETS"))
apoptosis_genes <- l1000.analysis::get_gene_set_genes(gsc = gsc, collection_name = "hallmark", gene_set_name = "HALLMARK_APOPTOSIS")

genes_of_interest <- intersect(colnames(merged_data_collapsed), c(cell_cycle_genes, apoptosis_genes))



n_fold_cv <- 10
folds <- cvTools::cvFolds(n = length(unique(merged_data_collapsed$pert_iname)), K = n_fold_cv)
unique_perts <- unique(merged_data_collapsed$pert_iname)
lambda_val <- 0.001
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

for (fold_ii in seq(n_fold_cv)) {
  
  ## split data
  curr_fold_test_idx <- which(folds$which == fold_ii)
  
  curr_train_perts <- unique_perts[-curr_fold_test_idx]
  curr_test_perts <- unique_perts[curr_fold_test_idx]
  
  curr_train_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, genes_of_interest] %>% as.matrix()
  curr_train_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, "auc_avg"] %>% as.matrix()
  
  curr_test_data <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, genes_of_interest] %>% as.matrix()
  curr_test_labels <- merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  lm_fit <- lm(formula = auc_avg ~ tas, data = merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_train_perts, c("tas", "auc_avg")])
  lm_preds <- predict(lm_fit, newdata = merged_data_collapsed[merged_data_collapsed$pert_iname %in% curr_test_perts, c("tas"), drop = F])
  curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
}

lm_pred_top_global_genes_df <- data.frame(
  lm_type = c(
    rep("expression_kfold", length(truth_pred_cor_list)),
    rep("tas_kfold", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double(),
    tas_results = truth_pred_cor_list_tas %>% as.double()
  )
)



combined_results <- rbind(
  lm_pred_df %>% mutate(gene_type = "lm"), 
  lm_pred_top_global_genes_df %>% mutate(gene_type = "top_global")
)
combined_results_melted <- melt(combined_results)
combined_results_melted <- combined_results_melted[grepl("expression", combined_results_melted$lm_type),]
combined_results_melted %>% ggplot(aes(gene_type, value)) + 
  geom_boxplot() + 
  cdsr::theme_Publication()







### COMPARE TO MODEL TRAINED ON AVERAGED EXPRESSION --------

# get average profile for each compound

mean_profile_list <- list()
auc_list <- list()
pert_iname_list <- list()
for (curr_cpd in unique(merged_data_collapsed$pert_iname)) {
  
  curr_profiles <- merged_data_collapsed[merged_data_collapsed$pert_iname == curr_cpd,]
  curr_mean_profile <- colMeans(curr_profiles[,lm_genes])
  
  mean_profile_list <- rbind(mean_profile_list, t(replicate(nrow(curr_profiles), curr_mean_profile)))
  
  
  auc_list <- c(auc_list, curr_profiles$auc_avg)
  pert_iname_list <- c(pert_iname_list, curr_profiles$pert_iname)
  
  
  
}

mean_profile_df <- mean_profile_list %>% as.data.frame() %>% 
  set_colnames(lm_genes)
mean_profile_df['pert_iname'] <- pert_iname_list %>% as.character()
mean_profile_df['auc_avg'] <- auc_list %>% as.double()


# train on mean-collapsed, test on uncollapsed


unique_perts <- unique(data_for_lm$pert_iname)
lambda_val <- 0.03
truth_pred_cor_list <- list()
truth_pred_cor_list_tas <- list()

pb = txtProgressBar(min = 0, max = length(unique_perts), initial = 0) 

for (cpd_ii in seq(length(unique_perts))) {
  
  
  ## split data
  curr_test_cpd <- unique_perts[cpd_ii]
  curr_train_cpds <- unique_perts[-which(unique_perts == curr_test_cpd)]
  
  # use full data for training
  curr_train_data <- mean_profile_df[mean_profile_df$pert_iname %in% curr_train_cpds, lm_genes] %>% as.matrix()
  curr_train_labels <- mean_profile_df[mean_profile_df$pert_iname %in% curr_train_cpds, "auc_avg"] %>% as.matrix()
  
  # use only selective compounds for testing
  curr_test_data <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, lm_genes] %>% as.matrix()
  curr_test_labels <- data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, "auc_avg"] %>% as.matrix()
  
  ## fit model
  glmnet_fit <- glmnet::glmnet(x = curr_train_data, y = curr_train_labels, lambda = lambda_val)
  
  ## get predictions
  curr_preds <- predict(glmnet_fit, newx = curr_test_data)
  
  ## check correlation of preds with truth
  curr_truth_pred_cor <- cor(curr_preds %>% as.double(), curr_test_labels %>% as.double())
  truth_pred_cor_list <- c(truth_pred_cor_list, curr_truth_pred_cor)
  
  
  
  ## FIT TAS MODEL
  # lm_fit <- lm(formula = auc_avg ~ tas, data = mean_profile_df[mean_profile_df$pert_iname %in% curr_train_cpds, c("tas", "auc_avg")])
  # lm_preds <- predict(lm_fit, newdata = data_for_lm[data_for_lm$pert_iname %in% curr_test_cpd, c("tas"), drop = F])
  # curr_truth_pred_cor_tas <- cor(lm_preds %>% as.double(), curr_test_labels %>% as.double())
  # truth_pred_cor_list_tas <- c(truth_pred_cor_list_tas, curr_truth_pred_cor_tas)
  
  setTxtProgressBar(pb, cpd_ii)
  
  
}

lm_pred_per_cpd_mean_profiles_df <- data.frame(
  lm_type = c(
    rep("Expression", length(truth_pred_cor_list))
    # rep("TAS", length(truth_pred_cor_list_tas))
  ),
  results = c(
    truth_pred_cor_list %>% as.double()
    # tas_results = truth_pred_cor_list_tas %>% as.double()
  ),
  cpd = c(rep(unique_perts, 2))
          # rep(unique_perts, 2))
)

lm_pred_per_cpd_df$cpd %<>% as.character()
lm_pred_per_cpd_mean_profiles_df$cpd %<>% as.character()
per_cpd_results_combined <- merge(lm_pred_per_cpd_df %>% filter(lm_type == "Expression"), 
                                  lm_pred_per_cpd_mean_profiles_df, 
                                  by = "cpd", 
                                  suffixes = c("_original", "_collapsed"))

p <- per_cpd_results_combined %>% 
  filter(lm_type_collapsed == "Expression") %>% 
  ggplot(aes(results_original, results_collapsed)) + 
  geom_point(size = 3) + 
  geom_abline() + 
  cdsr::theme_Publication() + 
  l1000.analysis::theme_scatterplot() + 
  labs(x = "Model Pearson corr.,\nOriginal data", y = "Model Pearson corr.,\nAveraged expression profiles")
p
ggsave(file.path(plot_dir, "mean_collapsed_predictive_comparison2.png"), p, height = 10, width = 12)





