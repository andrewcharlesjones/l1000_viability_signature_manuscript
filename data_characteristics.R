library(dplyr)
library(plyr)
library(magrittr)
library(taigr)
library(tibble)
library(limma)
library(ggplot2)

plot_dir <- "~/Documents/cancer_data_science/l1000_via_sig_ms/figures/plots/"


### LOAD DATA ----------

l1000_data <- l1000.analysis::load_l1000_cpd_data(normalize_sensitivity = T)
merged_data <- cbind(l1000_data$expression, l1000_data$sensitivity, l1000_data$metadata)
# 
# merged_data_collapsed <- merged_data %>% 
#   group_by(pert_iname, cell_id) %>% 
#   summarise_all(mean)

metadata_collapsed <- l1000_data$metadata %>% 
  group_by(pert_iname, cell_id, pert_dose, pert_time) %>% 
  summarise_all(mean)



### TAS BY TIME ------

metadata_collapsed %>% ggplot(aes(tas)) + 
  geom_density() + 
  cdsr::theme_Publication() + 
  facet_wrap( ~ pert_time) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")



### TAS BY DOSE ------

dose_counts <- table(metadata_collapsed$pert_dose)
dose_counts <- dose_counts[order(-dose_counts)]
common_doses <- dose_counts[1:9] %>% names()

metadata_collapsed %>% 
  filter(pert_dose %in% common_doses) %>% 
  ggplot(aes(tas)) +
  geom_density() + 
  cdsr::theme_Publication() + 
  facet_wrap( ~ pert_dose) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")


### PER-CPD VARIANCE OF TAS BY DOSE -----


dose_counts <- table(metadata_collapsed$pert_dose)
dose_counts <- dose_counts[order(-dose_counts)]
common_doses <- dose_counts[1:5] %>% names()

variance_table <- metadata_collapsed %>% 
  group_by(pert_iname, pert_dose, pert_time) %>% 
  summarise_all(var)

variance_table %>% 
  filter(pert_dose %in% common_doses) %>% ggplot(aes(tas)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity") + 
  geom_density() + 
  facet_grid(rows = vars(pert_time), cols = vars(pert_dose)) + 
  cdsr::theme_Publication()


### 2D DENSITY TAS BY SENSITIVITY -----

deduplicated_merged <- merged_data[!duplicated(merged_data[,c("cell_id", "pert_iname")]),]
deduplicated_merged %>% ggplot(aes(tas, auc_avg)) + 
  geom_point() + 
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
  cdsr::theme_Publication()


variance_table <- deduplicated_merged[,c("pert_iname", "pert_dose", "pert_time", "auc_avg", "cell_id", "tas")] %>% 
  group_by(pert_iname) %>% 
  summarise_all(var)

variance_table %>% ggplot(aes(tas, auc_avg)) + 
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
  cdsr::theme_Publication()
   
