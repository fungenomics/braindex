library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(dplyr)
library(feather)
library(glue)
library(ggplot2)
source("../functions.R")

# --------------------color palettes-------------------------

message("@ color palettes")

# make color palettes using the cluster label column and the colour label in the metadata data frame
# color palette for ggplot plots are named vectors with cluster colour as names and colour hex codes as values
# colour palette for pheatmap is a list of one element named "Cluster", this element is a named vector

#different metadata files are loaded here because different SCENIC datasets (non-extended, extended, per-time point) 
#have cluster labels generated at different dates 

#metadata <- read_tsv("shared/metadata_20190716.tsv") #for the non-extended dataset

metadata_extended <- read_tsv("shared/metadata_20210710_with_qc.tsv") #for the extended dataset in the joint space

#maps the per sample cluster label to a more broad cluster label
#F-e10_VRGC maps to RGC
#used in timeseries 
lvl2_cluster_extended <- metadata_extended %>% select(Label, Level2_type) %>% deframe

#even borader cluster mapping, for use in the heatmap
lvl1_cluster_extended <- metadata_extended %>% select(Label, Level1_type) %>% deframe

metadata_per_sample <- read_tsv("shared/metadata_20201028_with_qc.tsv") #for the per-timepoint analyses 

#used in timeseries 
lvl2_cluster_per_sample <- metadata_per_sample %>% select(Label, Level2_type) %>% deframe

#even borader cluster mapping, for use in the heatmap
lvl1_cluster_per_sample <- metadata_per_sample %>% select(Label, Level1_type) %>% deframe

cell_ontological_class_labels_from_lvl2 <- metadata_per_sample %>% filter(!is.na(Level2_type)) %>%
  select(Level2_type, Cell_ontological_class) %>% deframe

cell_ontological_class_labels_from_per_sample_cluster <- metadata_per_sample %>%
  select(Label, Cell_ontological_class) %>% deframe


lvl2_cluster_labels <- c(lvl2_cluster_extended, lvl2_cluster_per_sample)

lvl1_cluster_labels <- c(lvl1_cluster_extended, lvl1_cluster_per_sample)

cell_onto_label <- c(cell_ontological_class_labels_from_per_sample_cluster, 
                     cell_ontological_class_labels_from_lvl2)
#per_sample data colour palette
colour_palette_per_sample <- metadata_per_sample %>% select(Label, Colour) %>% deframe()
colour_palette_per_sample_space <- colour_palette_per_sample
names(colour_palette_per_sample_space) <- gsub("_", " ", names(colour_palette_per_sample))

# color palette for heatmap
colour_palette_cluster <- metadata_extended %>% 
  mutate(Label = gsub("_EXCLUDE", "", Label)) %>%
  # use gsub to change all contents in Cluster (cluster name format)
  mutate(Cluster = gsub("_", " ", Label)) %>%
  # Get two columns
  select(Cluster, Colour) %>%
  # Convert to vector of colours, where the first column gives the names
  # and the second column is converted to the values
  deframe() # VECTOR , not data frame

colour_palette_cluster_underscore <- colour_palette_cluster
names(colour_palette_cluster_underscore) <- gsub(" ", "_", names(colour_palette_cluster))


# color palette for timeseries plot, tab3
colour_palette <- metadata_extended %>% 
  mutate(Cluster = gsub("_", " ", Label)) %>% 
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

#------------------------tf_input-------------------
#this is list of the TF input into SCENIC, all active TF identified in each dataset is a subset of this
all_tf_list <- scan("shared/Mus_musculus_TF_one_TF_per_line.txt", character())


message("@ extended metadata")

# -----------------------joint_cortex_extended-----------------------------------
#joint_cortex_extended metadata containing DR coordinates 
forebrain_data_extended <- read_tsv("joint_cortex_extended/joint_cortex_extended.metadata.tsv", 
                                    col_types = cols(.default = col_character(), 
                                                     Joint_extended_PC_1 = col_double(),
                                                     Joint_extended_PC_2 = col_double(),
                                                     Joint_extended_tSNE_1 = col_double(),
                                                     Joint_extended_tSNE_2 = col_double(),
                                                     Joint_extended_UMAP_1 = col_double(),
                                                     Joint_extended_UMAP_2 = col_double())) %>% 
  transmute("Cell" = Cell_clean, "Sample" = Sample, "Sample_cluster" = ID_20210710_with_exclude,
            "Joint_cluster_number" = Joint_extended_cluster, "Joint_cluster" = ID_20210710_joint_clustering,
            "PC1" = Joint_extended_PC_1, "PC2" = Joint_extended_PC_2, 
            "tSNE_1" = Joint_extended_tSNE_1, "tSNE_2" = Joint_extended_tSNE_2,
            "UMAP1" = Joint_extended_UMAP_1, "UMAP2" = Joint_extended_UMAP_2) #renames different columns to be more clear

#TFs that are active in this dataset
TF_active_cortex_extended <- as_tibble(read_rds("joint_cortex_extended/joint_cortex_extended.active_regulons.Rds"))

#Dataframe containing TF-target gene pair in each row with info from the GRNBoost portion of SCENIC
#Use: GRN tab, Gene target info Table Tab
TF_target_gene_joint_cortex_extended <- as_tibble(read_rds(
  "joint_cortex_extended/joint_cortex_extended.regulon_target_info.Rds")) %>%
  select(-logo)

unique_TF_cortex_extended <- unique(TF_target_gene_joint_cortex_extended[["TF"]])

#Converts between TF's gene symbol and the name of its regulon produced by SCENIC
#regulon name format: {gene symbol}_{extended} ({number of genes in the regulon})
#extended regulon indicates if low confidence regulatory relationships are included in the regulon
# ex: Gene symbol: Arx   Regulon name: Arx_extended (21g)
TF_and_ext_cortex_extended <- identify_tf(TF_active_cortex_extended)

# binarized activity matrix for time-series ribbon plot
binary_activity_cortex_extended <- readRDS("joint_cortex_extended/4.2_binaryRegulonActivity_nonDupl_cortex_extended.Rds")
colnames(binary_activity_cortex_extended) <- gsub("_.","", colnames(binary_activity_cortex_extended))
tf_df_cortex_extended <- as_tibble(rownames(binary_activity_cortex_extended))

#metadata file for time-series ribbon plot
timeseries_input_meta_cortex_extended <- create_metadata_timeseries(forebrain_data_extended, "cortex", lvl2_cluster_labels)


data_cortex_extended <- list(
  "cell_metadata"  = forebrain_data_extended,
  "TF_and_ext" = TF_and_ext_cortex_extended,
  "TF_target_gene_info" = TF_target_gene_joint_cortex_extended,
  "unique_active_TFs_bare" = unique_TF_cortex_extended,
  "active_TFs" = TF_active_cortex_extended,
  "timeseries_input_meta" = timeseries_input_meta_cortex_extended,
  "binary_active_TFs" = tf_df_cortex_extended,
  "binary_activity" = binary_activity_cortex_extended
  
)

save(data_cortex_extended, file = "joint_cortex_extended/cortex_extended_prep.Rda")

# -----------------------joint_pons_extended-----------------------------------
#every element identical to the cortex_extended data
pons_data_extended <- read_tsv("joint_pons_extended/joint_pons_extended.metadata.tsv", 
                                    col_types = cols(.default = col_character(), 
                                                     Joint_extended_PC_1 = col_double(),
                                                     Joint_extended_PC_2 = col_double(),
                                                     Joint_extended_tSNE_1 = col_double(),
                                                     Joint_extended_tSNE_2 = col_double(),
                                                     Joint_extended_UMAP_1 = col_double(),
                                                     Joint_extended_UMAP_2 = col_double())) %>% 
  transmute("Cell" = Cell_clean, "Sample" = Sample, "Sample_cluster" = ID_20210710_with_exclude,
            "Joint_cluster_number" = Joint_extended_cluster, "Joint_cluster" = ID_20210710_joint_clustering,
            "PC1" = Joint_extended_PC_1, "PC2" = Joint_extended_PC_2, 
            "tSNE_1" = Joint_extended_tSNE_1, "tSNE_2" = Joint_extended_tSNE_2,
            "UMAP1" = Joint_extended_UMAP_1, "UMAP2" = Joint_extended_UMAP_2)

TF_active_pons_extended <- as_tibble(read_rds("joint_pons_extended/joint_pons_extended.active_regulons.Rds"))

TF_target_gene_joint_pons_extended <- as_tibble(read_rds(
  "joint_pons_extended/joint_pons_extended.regulon_target_info.Rds")) %>%
  select(-logo)

unique_TF_pons_extended <- unique(TF_target_gene_joint_pons_extended[["TF"]])

TF_and_ext_pons_extended <- identify_tf(TF_active_pons_extended)

timeseries_input_meta_pons_extended <- create_metadata_timeseries(pons_data_extended, "pons", lvl2_cluster_labels)

binary_activity_pons_extended <- readRDS("joint_pons_extended/4.2_binaryRegulonActivity_nonDupl_pons_extended.Rds")
colnames(binary_activity_pons_extended) <- gsub("_.","", colnames(binary_activity_pons_extended))
tf_df_pons_extended <- as_tibble(rownames(binary_activity_pons_extended))



data_pons_extended <- list(
  "cell_metadata"  = pons_data_extended,
  "TF_and_ext" = TF_and_ext_pons_extended,
  "TF_target_gene_info" = TF_target_gene_joint_pons_extended,
  "unique_active_TFs_bare" = unique_TF_pons_extended,
  "active_TFs" = TF_active_pons_extended,
  "timeseries_input_meta" = timeseries_input_meta_pons_extended,
  "binary_active_TFs" = tf_df_pons_extended,
  "binary_activity" = binary_activity_pons_extended
  
)
save(data_pons_extended, file = "joint_pons_extended/pons_extended_prep.Rda")

message("@ master color palettes")

#------------------------master color palette----------------------------------
#When I was adding various different data sets generated at different times with different formats,
#colour palettes were super annoying. 
#The master palette contains colours for all cluster labels regardless of when the cluster labels were assigned or
#if the labels have underscores or not
#Used in all ggplots in the app

#Note: all cluster labels should have the same format now, but will keep using this in case I missed something

extended_mouse_joint_cluster_palette <- readRDS("shared/palette_ID_20210710_joint_clustering.Rds")

#palette for the level2 cluster labels from metadata_extended
palette_broad_clusters <- c("RGC" = "#ffcc00",
                            "Telencephalic progenitors" = "#FFE5AF",
                            "Inhibitory neurons" = "#135ca0",
                            "Neuronal IPC" = "#f77750",
                            "Excitatory neurons" = "#840200",
                            "Meninges" = "#dbd2d7",
                            "Unresolved" = "gray50", 
                            "Myeloid" = "gray50",
                            "Endothelial" = "#636363",
                            "Thalamic precursors" = "gray50",
                            "Other neurons" = "#c18ba0",
                            "Cortical hem" = "#FFE5AF",
                            "Thalamic neurons" = "#805f91",
                            "Pericytes" = "#857e89",
                            "Microglia" = "#aca2b2",
                            "Gliogenic progenitors" = "#d5d98b",
                            "Astrocytes" = "#00a385",
                            "Choroid plexus" = "#1ba02a",
                            "Oligodendrocytes" = "#d5d98b",
                            "Macrophages" = "#86778e",
                            "Ependymal" = "#8ee5cf",
                            "Neurons" = "#ff9385",
                            "Schwann cells" = "#4b6d34",
                            "Vascular smooth muscle" = "#665f5f",
                            "Hindbrain progenitors" = "#ffe5af",
                            "Vascular leptomeningeal" = "#ceb9c5",
                            "Mixed progenitors" = "#ffe500"
                            
)

lvl1_cluster_palette <- c("Progenitors" = "#ffbda3",
                                  "Neurons"    = "#135ca0",
                                  "Leptomeningeal" = "#ceb9c5",
                                  "Unresolved" = "gray50",
                                  "Blood" = "gray90",
                                  "Endothelial" = "#636363",
                                  "Vascular" = "#665f5f",
                                  "Immune" = "#86778e",
                                  "Glia" = "#00a385")

cell_ontological_class_palette <- c("RGC"                  = "#ffcc00",
                                    "Glial progenitors"    = "#d5d98b",
                                    "OPC"                  = "#e0de53",
                                    "Proliferating OPC"    = "#e6f957",
                                    "Oligodendrocytes"     = "#b4e04e",
                                    "Astrocytes"           = "#00a385",
                                    "Ependymal"            = "#8ee5cf",
                                    "Neuronal progenitors" = "#ffbda3",
                                    "Neurons"              = "#135ca0",
                                    "Immune"               = "gray50",
                                    "Vascular & other"     = "gray70",
                                    "Normal"               = "gray90")


master_palette <- c(hm_anno_new$side_colors$Cluster,        # per-timepoint cluster with timepoint removed
                    colour_palette_cluster,                 # per-timepoint cluster, with spaces
		                colour_palette_cluster_underscore,      # per-timepoint cluster, with underscores
		               #forebrain_cluster_palette,              # joint clustering, forebrain
		               #pons_cluster_palette,                   # joint clustering, pons
                    extended_mouse_joint_cluster_palette,
		                colour_palette_per_sample_space,
		                colour_palette_per_sample)

master_palette <- list("Cluster" = master_palette, 
                       "Broad Cluster" = palette_broad_clusters,
                       "Broader Cluster" = lvl1_cluster_palette,
                       "Cell Ontology Class" = cell_ontological_class_palette)

message("@ time point data")

#---------------------time_point data----------------------------------------------
#Data processing for each time point is the same
for (reg in c("ct", "po")){
  
  for(tp in c("e10", "e12", "e13", "e15", "e16", "e18", "p0", "p3", "p6")){
    TF_target_gene_info <- as_tibble(read_rds(glue("{reg}_{tp}/{reg}_{tp}.regulon_target_info.Rds"))) %>%
      select(-logo)
    
    unique_TF <- unique(TF_target_gene_info[["TF"]])
    
    TF_active <- as_tibble(read_rds(glue("{reg}_{tp}/{reg}_{tp}.active_regulons.Rds")))
    TF_and_ext <- identify_tf(TF_active)
    
    #metadata with DR coords 
    cell_data <- read_tsv(glue("{reg}_{tp}/{reg}_{tp}.metadata.tsv")) 
    
    #indicates cells that should not be included in plots because they belong to a blacklisted cluster
    #Used to filter the TF activity per cell feather in the create_activity_data function 
    black_list_cells <- cell_data %>% select(cell, ID_20201028_with_exclude) %>%
      filter(grepl("EXCLUDE", ID_20201028_with_exclude)) %>% select(cell) %>%
      deframe()
    
    cell_data <- cell_data %>% filter(!grepl("EXCLUDE", ID_20201028_with_exclude))
      
    
    x <- list(
           "TF_target_gene_info" = TF_target_gene_info,
           "TF_and_ext" = TF_and_ext,
           "cell_metadata" = cell_data,
           "bad_cells" = black_list_cells,
           "unique_TF" = unique_TF
    )
    
    saveRDS(x, file = glue("{reg}_{tp}/{reg}_{tp}_prep.Rds"))
  }
  
}

message("@ shared")

# -----------------------------shared data-----------------------------
save(colour_palette_cluster,
     hm_anno, hm_anno_new, colour_palette, all_tf_list, 
     master_palette, palette_broad_clusters, 
     lvl1_cluster_labels, lvl2_cluster_labels, cell_onto_label,
     file = "shared/common_prep.Rda")

message("@ ribbon")

#-----------------cell proportion over time ribbon plot--------------------

#calculates the fraction of cells of the total that belong to each cluster
forebrain_fraction <- forebrain_data_extended %>% select(Cell, Sample, Sample_cluster) %>%
  filter(!grepl("EXCLUDE", Sample_cluster))

forebrain_fraction <- forebrain_fraction %>% 
  mutate(broad_cluster = recode(forebrain_fraction$Sample_cluster, !!!lvl2_cluster_labels)) %>% 
  separate(Sample_cluster, into = c("tp", "Cluster"), sep = "_") %>%
  group_by(Sample) %>% mutate (total_in_tp = n()) %>%
  ungroup() %>% group_by(Sample, broad_cluster) %>%
  mutate(frac = n()/total_in_tp) %>% 
  ungroup() %>% group_by(Sample)


forebrain_clusters <- forebrain_fraction$broad_cluster %>% unique

#removes duplicate rows such that each row corresponds to data for one cluster at one timepoint
unique_forebrain_fraction <- forebrain_fraction %>% select(-Cell, -Cluster) %>% 
  distinct() %>% select(-total_in_tp, -tp)


#not all clusters appear in all timepoints, each timepoint needs to have a complete set of clusters for ribbon plot to look right
#adds any missing clusters in each timepoint and sets its fraction value to 0
unique_forebrain_fraction_complete <- unique_forebrain_fraction %>%
  mutate(broad_cluster = factor(broad_cluster, levels = unique(.$broad_cluster))) %>%
  complete(broad_cluster, nesting(Sample), fill = list(frac = 0))


unique_forebrain_fraction_complete$xpos = group_indices(unique_forebrain_fraction_complete)


forebrain_plot <- unique_forebrain_fraction_complete %>%
  ggplot(aes(x = xpos, y = frac, fill = broad_cluster)) +
  geom_area(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = palette_broad_clusters, drop = FALSE, name = "") +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                     labels = c("E10.5", "E12.5", "E13.5", "E15.5", "E16.5", "E18.5", "P0", "P3", "P6"),
                     limits = c(1, 9)) +
  labs(x = "Developmental Age", y = "Proportion") +
  guides(fill = guide_legend(ncol = 5)) +
  theme_min() + 
  theme(legend.position = "bottom")

#same for pons


pons_fraction <- pons_data_extended %>% select(Cell, Sample, Sample_cluster) %>%
  filter(!grepl("EXCLUDE", Sample_cluster))

# for pons metadata, some cells are labelled E15.5 for sample, and E12.5 for the cluster labels
# making the assumption here that the sample label is right and the cluster label is wrong
# in the code: group_by(Sample) and doing the total cell and proportion calculations based on assumption that 
# the sample label is correct. All these confused cells are put with the E15.5 cells.
pons_fraction <- pons_fraction %>%
  mutate(broad_cluster = recode(pons_fraction$Sample_cluster, !!!lvl2_cluster_labels)) %>%
  separate(Sample_cluster, into = c("tp", "Cluster"), sep = "_") %>%
  group_by(Sample) %>% mutate (total_in_tp = n()) %>% 
  ungroup() %>% group_by(Sample, broad_cluster) %>%  
  mutate(frac = n()/total_in_tp, number = n()) %>% 
  ungroup() %>% group_by(Sample)


unique_pons_fraction <- pons_fraction %>% select(-Cell, -Cluster) %>% 
  distinct() %>% select(-total_in_tp, -tp) 


unique_pons_fraction_complete <- unique_pons_fraction %>%
  mutate(broad_cluster = factor(broad_cluster, levels = unique(.$broad_cluster))) %>%
  complete(broad_cluster, nesting(Sample), fill = list(frac = 0))


unique_pons_fraction_complete$xpos = group_indices(unique_pons_fraction_complete)


pons_plot <- unique_pons_fraction_complete %>%
  ggplot(aes(x = xpos, y = frac, fill = broad_cluster)) +
  geom_area(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = palette_broad_clusters, drop = FALSE, name = "") +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                     labels = c("E10.5", "E12.5", "E13.5", "E15.5", "E16.5", "E18.5", "P0", "P3", "P6"),
                     limits = c(1, 9)) +
  labs(x = "Developmental Age", y = "Proportion") +
  guides(fill = guide_legend(ncol = 5)) +
  theme_min() + 
  theme(legend.position = "bottom")

timeseries_proportion_plots <- list(
  "forebrain_plot" = forebrain_plot,
  "pons_plot" = pons_plot
)

save(timeseries_proportion_plots, file = "shared/timeseries_proportion_plots.Rda")

message("@ done.")

