library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(dplyr)
library(feather)
library(glue)
library(ggplot2)
source("../functions.R")

# ———————————————————————————————————color palette————————————————————————————————————————
# make color palette
metadata <- read_tsv("shared/metadata_20190716.tsv")

# color palette for heatmap
colour_palette_cluster <- metadata %>%
  # use gsub to change all contents in Cluster (cluster name format)
  mutate(Cluster = gsub("_", " ", Cluster)) %>%
  # Get two columns
  select(Cluster, Colour) %>%
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe() # VECTOR , not data frame

all_tf_list <- scan("shared/Mus_musculus_TF_one_TF_per_line.txt", character())

# color palette for timeseries plot, tab3
colour_palette <- metadata %>% 
  mutate(Cluster = gsub("_", " ", Cluster)) %>% 
  separate(Cluster, into = c("Prefix", "Cluster"), sep = " ") %>% 
  # Get two columns
  select(Cluster, Colour) %>% 
  distinct(Cluster, .keep_all = TRUE) %>% 
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe()

# A helper function to prepare a dataframe to annotate the heatmap with colours
hm_anno <- makePheatmapAnno(colour_palette_cluster, "Cluster")
# this is used for generating the anno_row, since we need to have the same row names as
# those in t(act_cluster) to color correctly

hm_anno_new <- makePheatmapAnno(colour_palette, "Cluster")
# this is used in: annotation_colors = hm_anno_new$side_colors, in both heatmaps (by cluster/cells)



# ———————————————————————————————————Cortex data————————————————————————————————————————
forebrain_data <- read_tsv("joint_cortex/Forebrain_join.2D.tsv") %>% # for UMAP cluster
  mutate(Sample_cluster = str_replace(Sample_cluster," ","_"))
# clean some samples with space in between ...

TF_active <- as_tibble(read_rds("joint_cortex/joint_cortex.active_regulons.Rds"))

# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene <- as_tibble(read_rds("joint_cortex/joint_cortex.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF <- unique(TF_target_gene[["TF"]])

#reads metadata file for color palette of clustering by region
forebrain_cluster_palette <- read_tsv("joint_cortex/Jessa2019_Table_2b_joint_cortex_metadata.tsv")
forebrain_cluster_palette <- forebrain_cluster_palette %>% select(Cluster, Colour) %>% deframe()

TF_and_ext <- identify_tf(TF_active)

timeseries_input_meta_cortex <- create_metadata_timeseries(forebrain_data, "cortex")


# metadata specific for each cell, corresponding to the activity data
#cell_metadata_cortex_prep <- read_tsv("joint_cortex/joint_cortex.metadata.tsv")

#cell_metadata_cortex_test <- create_cell_metadata_cortex(forebrain_data)

#cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex_prep)
# activity for cortex timeseries graph data
binary_activity <- readRDS("joint_cortex/joint_cortex.binaryRegulonActivity_nonDupl.Rds")
tf_df <- as_tibble(rownames(binary_activity)) #a dataframe that contains all the tf with 
# best representation of its identity in the binary_activity dataset

l <- c()
l_nexist_cortex <- c()
for (tf in unique_TF){
  tf_after <- translate_tf(tf, tf_df)
  if(tf_after !=FALSE ){
    l <- c(l, tf)
  }
  else{l_nexist_cortex<- c(l_nexist_cortex,tf)}
}

# ----------------------------------Pons data-------------------------------------------------------

pons_data <- read_tsv("joint_pons/Pons_join.2D.tsv") # for UMAP cluster

TF_active_pon <- as_tibble(read_rds("joint_pons/joint_pons.active_regulons.Rds"))

# These datasets describe TF and genes that are target of TFs, don't have ext suffix
TF_target_gene_pon <- as_tibble(read_rds("joint_pons/joint_pons.regulon_target_info.Rds")) %>%
  select(-logo)
unique_TF_pon <- unique(TF_target_gene_pon[["TF"]])

#reads metadata file for color palette of clustering by region
pons_cluster_palette <- read_tsv("joint_pons/Jessa2019_Table_2c_joint_pons_metadata.tsv")
pons_cluster_palette <- pons_cluster_palette %>% select(Cluster, Colour) %>% deframe()


TF_and_ext_pon <- identify_tf(TF_active_pon)

timeseries_input_meta_pons <- create_metadata_timeseries(pons_data,"pons") %>%
  filter(Cell != "___po_e12_TACGGGCGTCAAGCGA")
# filter out the extra cell
# remove the extra line to make the number of cells the same as the binary activity pon data
# to correctly make the timeseires ribbon plot

# activity for cortex timeseries graph data
binary_activity_pon <- readRDS("joint_pons/joint_pons.binaryRegulonActivity_nonDupl.Rds")
tf_df_pon <- as_tibble(rownames(binary_activity_pon)) #a dataframe that contains all the tf with 
# best representation of its identity in the binary_activity dataset

l <- c()
l_nexist_pons <- c()
for (tf in unique_TF_pon){
  tf_after <- translate_tf(tf, tf_df_pon)
  if(tf_after !=FALSE ){
    l <- c(l, tf)
  }
  else{l_nexist_pons<- c(l_nexist_pons,tf)}
}


# make two lists containing same name (will be assigned to a reactive list),
# then we can use the same name to code
data_cortex <- list(
  "cell_metadata"  = forebrain_data,
  "TF_and_ext" = TF_and_ext,
  "TF_target_gene_info" = TF_target_gene,
  "unique_active_TFs_bare" = unique_TF,
  "active_TFs" = TF_active,
  "binary_active_TFs" = tf_df,
  "timeseries_input_meta" = timeseries_input_meta_cortex,
  "binary_activity" = binary_activity,
  "tfs_not_exist_timeseries" = l_nexist_cortex,
  "cluster_palette" = forebrain_cluster_palette

)

data_pons <- list(
  "cell_metadata" = pons_data,
  "TF_and_ext" = TF_and_ext_pon,
  "TF_target_gene_info" = TF_target_gene_pon,
  "unique_active_TFs_bare" = unique_TF_pon,
  "active_TFs" = TF_active_pon,
  "binary_active_TFs" = tf_df_pon,
  "timeseries_input_meta" = timeseries_input_meta_pons,
  "binary_activity" = binary_activity_pon,
  "tfs_not_exist_timeseries" = l_nexist_pons,
  "cluster_palette" = pons_cluster_palette
)

master_palette <- append(hm_anno_new$side_colors$Cluster, forebrain_cluster_palette)
master_palette <- append(master_palette, pons_cluster_palette)
master_palette <- list("Cluster" = master_palette)
#---------------------ct_e12 data----------------------------------------------
#use a loop for this 
for (reg in c("ct", "po")){
  
  for(tp in c("e12", "e15", "p0", "p3", "p6")){
    TF_target_gene_info <- as_tibble(read_rds(glue("{reg}_{tp}/{reg}_{tp}.regulon_target_info.Rds"))) %>%
      select(-logo)
    
    unique_TF <- unique(TF_target_gene_info[["TF"]])
    
    TF_active <- as_tibble(read_rds(glue("{reg}_{tp}/{reg}_{tp}.active_regulons.Rds")))
    TF_and_ext <- identify_tf(TF_active)
    
    cell_data <- read_tsv(glue("{reg}_{tp}/{reg}_{tp}.metadata.tsv")) 
    
    black_list_cells <- cell_data %>% select(Cell, ID_20190715_with_blacklist) %>%
      filter(grepl("BLACKLISTED", ID_20190715_with_blacklist)) %>% select(Cell) %>%
      deframe()
    
    #print(black_list_cells)
    
    cell_data <- cell_data %>% filter(!grepl("BLACKLISTED", ID_20190715_with_blacklist))
      
    
    awoo <- switch(reg, "ct" = "F-", "po" = "P-")
    
    dr_palette <- metadata %>% 
      separate(Cluster, into = c("Timepoint", "Cluster"), sep = "_") %>%
      filter(Timepoint == glue("{awoo}{tp}")) %>%
      unite(col = "Cluster", c("Timepoint", "Cluster"), sep = "_") %>%
      select(Cluster, Colour) %>% 
      deframe()
    
    x <- list(
           "TF_target_gene_info" = TF_target_gene_info,
           "TF_and_ext" = TF_and_ext,
           "cell_metadata" = cell_data,
           "bad_cells" = black_list_cells,
           "cluster_palette" = dr_palette,
           "unique_TF" = unique_TF
    )
    #print(x)
    
    saveRDS(x, file = glue("{reg}_{tp}/{reg}_{tp}_prep.Rds"))
  }
  
}

#-----------------------------forebrain joint cluster regulon activity data for heatmap------------
forebrain_regulon_activity_data <-
  read_feather("joint_cortex/joint_cortex.regulon_activity_per_cell.feather")

forebrain_joint_cluster_info <- forebrain_data %>% select(Cell, Joint_cluster)

forebrain_cluster_regulon_data <- 
  inner_join(forebrain_regulon_activity_data, forebrain_joint_cluster_info, by = "Cell")

forebrain_cluster_regulon_data <- forebrain_cluster_regulon_data %>% group_by(Joint_cluster) %>% 
  summarize_if(is.numeric, mean)

write_feather(forebrain_cluster_regulon_data, 
     path = "joint_cortex/joint_cortex.regulon_activity_per_joint_cluster.feather")

#-----------------------------pons joint cluster regulon activity data for heatmap------------
pons_regulon_activity_data <-
  read_feather("joint_pons/joint_pons.regulon_activity_per_cell.feather")

pons_joint_cluster_info <- pons_data %>% select(Cell, Joint_cluster)

pons_cluster_regulon_data <- 
  inner_join(pons_regulon_activity_data, pons_joint_cluster_info, by = "Cell")

pons_cluster_regulon_data <- pons_cluster_regulon_data %>% group_by(Joint_cluster) %>% 
  summarize_if(is.numeric, mean)

write_feather(pons_cluster_regulon_data, 
              path = "joint_pons/joint_pons.regulon_activity_per_joint_cluster.feather")
# ---------------------------cortex data-----------------------------
save(data_cortex, file = "joint_cortex/cortex_prep.Rda")


# -----------------------------pons data-----------------------------
save(data_pons, file = "joint_pons/pons_prep.Rda")


# -----------------------------shared data-----------------------------
save(colour_palette_cluster,
     hm_anno, hm_anno_new, colour_palette, all_tf_list, master_palette, file = "shared/common_prep.Rda")

#-----------------cell proportion over time ribbon plot--------------------

forebrain_fraction <- forebrain_data %>% select(Cell, Sample, Sample_cluster) %>%
  filter(!grepl("BLACKLIST", Sample_cluster)) %>%
  separate(Sample_cluster, into = c("tp", "Cluster"), sep = "_") %>%
  group_by(Sample) %>% mutate (total_in_tp = n()) %>%
  ungroup() %>% group_by(Sample, Cluster) %>%
  mutate(frac = n()/total_in_tp) %>% 
  ungroup() %>% group_by(tp)

forebrain_clusters <- forebrain_fraction$Cluster %>% unique

unique_forebrain_fraction <- forebrain_fraction %>% select(-Cell) %>% 
  unique %>% select(-Sample, -total_in_tp)

tp <- unique_forebrain_fraction$tp %>% unique()

for (i in tp){
  tp_clusters <- unique_forebrain_fraction %>% ungroup() %>%
    filter(tp == i) %>% pull(Cluster)
  #print(tp_clusters)
  clust_not_in <- forebrain_clusters[!(forebrain_clusters %in% tp_clusters)]
  #print(clust_not_in)
  to_add <- tibble(tp = rep(i, length(clust_not_in)), Cluster = clust_not_in, frac = rep(0, length(clust_not_in)))
  # for(j in clust_not_in){
  #   row <- tibble(tp = i, Cluster = j, frac = 0)
  # }
  unique_forebrain_fraction <- rbind(unique_forebrain_fraction, to_add)
}

unique_forebrain_fraction %>% ungroup() %>% group_by(tp)
unique_forebrain_fraction$xpos = group_indices(unique_forebrain_fraction)


forebrain_plot <- unique_forebrain_fraction %>%
  ggplot(aes(x = xpos, y = frac, fill = Cluster)) +
  geom_area(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = colour_palette, drop = FALSE, name = "") +
  scale_x_continuous(breaks = c(1,2,3,4,5),
                     labels = c("E12.5", "E15.5", "P0", "P3", "P6"),
                     limits = c(1, 5)) +
  labs(x = "Developmental Age", y = "Proportion") +
  guides(fill = guide_legend(ncol = 5)) +
  theme_min() + 
  theme(legend.position = "bottom")

#same for pons

pons_fraction <- pons_data %>% select(Cell, Sample, Sample_cluster) %>%
  filter(!grepl("BLACKLIST", Sample_cluster)) %>%
  separate(Sample_cluster, into = c("tp", "Cluster"), sep = "_") %>%
  group_by(tp) %>% mutate (total_in_tp = n()) %>%
  ungroup() %>% group_by(Sample, Cluster) %>%
  mutate(frac = n()/total_in_tp, number = n()) %>% 
  ungroup() %>% group_by(tp)

pons_clusters <- pons_fraction$Cluster %>% unique

unique_pons_fraction <- pons_fraction %>% select(-Cell) %>% 
  unique %>% select(-Sample, -total_in_tp) 

test <- unique_pons_fraction %>% group_by(tp) %>% summarize(sum(frac))

tp <- unique_pons_fraction$tp %>% unique()

for (i in tp){
  tp_clusters <- unique_pons_fraction %>% ungroup() %>%
    filter(tp == i) %>% pull(Cluster)
  #print(tp_clusters)
  clust_not_in <- pons_clusters[!(pons_clusters %in% tp_clusters)]
  #print(clust_not_in)
  to_add <- tibble(tp = rep(i, length(clust_not_in)), Cluster = clust_not_in, frac = rep(0, length(clust_not_in)))
  # for(j in clust_not_in){
  #   row <- tibble(tp = i, Cluster = j, frac = 0)
  # }
  unique_pons_fraction <- rbind(unique_pons_fraction, to_add)
}

unique_pons_fraction %>% ungroup() %>% group_by(tp)
unique_pons_fraction$xpos = group_indices(unique_pons_fraction)


pons_plot <- unique_pons_fraction %>%
  ggplot(aes(x = xpos, y = frac, fill = Cluster)) +
  geom_area(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = colour_palette, drop = FALSE, name = "") +
  scale_x_continuous(breaks = c(1,2,3,4,5),
                     labels = c("E12.5", "E15.5", "P0", "P3", "P6"),
                     limits = c(1, 5)) +
  labs(x = "Developmental Age", y = "Proportion") +
  guides(fill = guide_legend(ncol = 5)) +
  theme_min() + 
  theme(legend.position = "bottom")

timeseries_proportion_plots <- list(
  "forebrain_plot" = forebrain_plot,
  "pons_plot" = pons_plot
)

save(timeseries_proportion_plots, file = "shared/timeseries_proportion_plots.Rda")



