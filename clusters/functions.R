
# Load required packages ----
library(feather)
library(tidyr)
library(dplyr)
library(cowplot)
library(glue)
library(stringr)
library(ggplot2)

# Set-up / load common data ----

# Cluster-level metadata
metadata <- data.table::fread("data/joint_mouse/metadata_20190715_select.tsv",
                              data.table = FALSE)

# Red-blue gradient and brain region colour palettes
load("data/joint_mouse/palettes.Rda")

# Joint mouse colour palette
load("data/joint_mouse/joint_mouse.palette_ID_20190715.Rda")

# Vector specifying the order of clusters in the dendrogram
load("data/joint_mouse/ID_20190715_dendrogram_order.Rda")


# Custom functions ----

#' Prepare input for bubble_plot
#' 
#' Load gene expression data from feather, tidy & optionally scale expression,
#' and return the dataframe required as input for the bubble_plot() function
#' 
#' @param gene Character vector, one or more genes of interest to plot
#' @param scale Logical, whether or not to linearly scale gene expression across
#' clusters to [0,1] to improve visualization. Default: TRUE
#' TODO: Use this to provide an option to donwload the underlying data.
#' 
#' @examples 
#' bubble_prep("Dlx1")
bubble_prep <- function(gene,
                        scale = TRUE) {
  
  # Load the mean expression of genes across clusters, given gene of interest
  exp <- read_feather("data/joint_mouse/mean_expression_per_ID_20190715_cluster.feather",
                      columns = c("Cluster", gene)) %>% 
    filter(Cluster %in% dendrogram_order)
  
  # Scale expression of each gene linearly across clusters to [0, 1]
  if (scale) {
    
    exp <- exp %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = "Cluster") %>%
      apply(2, scales::rescale, to = c(0, 1)) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "Cluster")
    
  }
  
  # Convert to long / tidy format with columns: Cluster, Gene, Expression
  exp <- exp %>%
    gather(., "Gene", "Expression", 2:ncol(.))
  
  # Load the prorportion of cells in each cluster in which each gene was detected,
  # and convert to long / tidy format with columns: Cluster, Gene, Pct1
  pct1 <- read_feather("data/joint_mouse/pct1_per_ID_20190715_cluster.feather",
                       columns = c("Cluster", gene)) %>%
    gather(., "Gene", "Pct1", 2:ncol(.))
  
  # Join with cluster metadata
  df <- left_join(exp, pct1, by = c("Cluster", "Gene"))  %>% 
    left_join(metadata, by = c("Cluster" = "Cluster_nounderscore"))
  
  # Tidy data for plotting
  df <- df %>%
    
    # Order genes to match order input by user
    mutate(Gene = factor(Gene, levels = rev(gene))) %>% 
    arrange(Gene) %>% 
    
    # Pad gene names so that the plot takes up a more standardized
    # width; to roughly the the # of characters in the gene w/ longest name
    # However, letters take up more pixels thn spaces, so do less padding
    # for genes with longer names
    mutate(Gene_padded = case_when(
      str_length(Gene) <= 5 ~ str_pad(Gene, 15, side = 'right', pad = " "),
      str_length(Gene) > 5 ~ str_pad(Gene, 12, side = 'right', pad = " "))
    ) %>% 
    mutate(Gene_padded = factor(Gene_padded, levels = unique(.$Gene_padded))) %>% 
    
    # Order the clusters on the x-axis to match the dendrogram image
    mutate(Cluster = factor(Cluster, levels = dendrogram_order)) %>%
    
    filter(!is.na(Cluster)) %>% 
    
    # Convert NAs (undetected genes) to 0s -- this ensures all
    # clusters have a value for all genes, so that all clusters are plot,
    # even if the gene was undetected
    replace_na(list(Expression = 0, Pct1 = 0)) 
  
  return(df)
  
}


#' Bubbleplot of gene expression
#' 
#' Generate a bubble plot for genes of interest across clusters in the mouse
#' dendrogram, where bubble colour encodes the mean expression in each cluster
#' and bubble size encodes the proportion of cells where each gene is detected
#'
#' @param df Data frame as returned by bubble_prep(), with require columns Cluster,
#' Gene_padded, Pct1, and Expression
#'
#' @return ggplot2 object
#'
#' @examples
#' bubble_prep("Dlx1") %>% bubbleplot()
bubble_plot <- function(df, max_point_size) {
  
  # Generate plot
  p1 <- df %>% 
    ggplot(aes(x = Cluster, y = Gene_padded)) +
    geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
    scale_size_area(max_size = max_point_size) +
    scale_color_gradientn(colours = tail(rdbu, 70)) +
    theme_min() +
    ylab(NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     colour = joint_mouse_palette, size = rel(0.7)),
          panel.grid.major.x = element_line(colour = "grey90"),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          # Do not show the legend because it is included in the static
          # dendrogram image displayed above the bubbleplot
          legend.position = "bottom") +
    # Put gene labels on the right hand side to improve alignment
    scale_y_discrete(position = "right")
  
  return(p1)
  
}


#' Prepare input for ribbon plot
#' 
#' Loads gene expression values, cell metadata, and tidies values while filtering
#' out blacklisted clusters
#'
#' @param gene String corresponding to gene name (only a single value allowed) 
#' @param region String, corresponding to the brain region to plot, one of
#' "joint_cortex" or "joint_pons"
#'
#' @return Dataframe with first column "Cell" and the rest corresponding to 
#' sample & cluster-level metadata for each Cell
#'
#' TODO: At some point, generalize this to be able to be used by other functions
#' which require a similar input
#'
#' @examples
#' prep_ribbon_input("Pdgfra", "joint_cortex")
prep_ribbon_input <- function(gene, region) {
  
  ribbon_df <- read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                            c("Cell", "ID_20190715_with_blacklist_and_refined", gene)) %>% 
    rowwise() %>% 
    mutate(Cell = str_extract(Cell, "[ACGT]{16}")) %>% 
    select(Cell, Cluster = ID_20190715_with_blacklist_and_refined, everything()) %>% 
    ungroup() %>% 
    left_join(metadata, by = "Cluster") %>% 
    mutate(Age = ifelse(Age == "E12.5-E15.5", "E12.5", Age))
  
  ribbon_df <- ribbon_df %>% 
    # Filter out cell types to exclude
    filter(!grepl("BLACKLISTED", Cluster)) %>% 
    filter(!is.na(Cell_type)) %>% 
    mutate(Age = factor(Age, levels = c("E12.5",
                                        "E15.5",
                                        "P0",
                                        "P3",
                                        "P6"))) %>% 
    arrange(Age)
  
}



#' Generate a ribbon plot
#' 
#' This function generates a ribbon plot representing proportions over time,
#' as seen in ,. This function is based on R code provided by the authors at
#' https://github.com/MarioniLab/EmbryoTimecourse2018/tree/master/analysis_scripts/atlas/vis/ribbon
#'
#' @param gene String corresponding to gene name (only a single value allowed) 
#' @param region String, corresponding to the brain region to plot, one of
#' "joint_cortex" or "joint_pons"
#' @param ymax Numeric, value in [0, 1] specifying the maximum value for the y-axis.
#' By default, y-axis is scaled to the range of the data (see more at 
#' https://ggplot2.tidyverse.org/reference/lims.html)
#'
#' @return ggplot2 object
#'
#' @examples
#' ribbon_plot("Pdgfra", "joint_cortex")
ribbon_plot <- function(gene, region, ymax = NA) {
  
  # Adapt palette to brain region
  if (region == "joint_cortex") colours <- cortex_palette
  else if (region == "joint_pons") colours <- pons_palette
  
  # Prep input dataframe
  ribbon_df <- prep_ribbon_input(gene, region)
  ribbon_df$gene <- ribbon_df[[gene]]
  
  # For each cluster at each timepoint, calculate the proportion of cells in
  # which the gene is detected
  ribbon_df_celltype_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    mutate(total = n()) %>% 
    group_by(Age, Cell_type) %>% 
    mutate(frac = sum(gene > 0) / total) %>% 
    distinct(Age, Cell_type, frac) %>% 
    ungroup()
  
  # For each timepoint, calculate the proportion of cells in which the gene 
  # is detected
  ribbon_df_cum_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    summarize(cumfrac = sum(gene > 0) / n()) %>% 
    ungroup()
  
  # Get the values to use on the x-axis
  timepoints2 <- ribbon_df$Age
  clusters <- ribbon_df$Cell_type
  
  df = data.frame(cluster = rep(unique(clusters), length(unique(timepoints2))),
                  stage = do.call(c, lapply(as.character(unique(timepoints2)), rep, times = length(unique(clusters)))))
  
  # Use the same order for clusters (vertically) as they are saved
  # in the colour palette
  df$ranking = match(df$cluster, names(colours))
  df = df[order(df$stage, df$ranking),]
  
  df <- left_join(df, select(ribbon_df_celltype_frac, cluster = Cell_type, stage = Age, frac)) %>% 
    # Complete cases when genes were not detected in certain timepoints/clusters
    # by replacing with a zero
    mutate(frac = replace_na(frac, 0)) %>% 
    left_join(select(ribbon_df_cum_frac, stage = Age, cumfrac))
  
  df$xpos = match(df$stage, unique(timepoints2))
  
  p1 <- df %>%
    ggplot(aes(x = xpos, y = frac, fill = cluster)) +
    geom_area(stat = "identity") +
    scale_fill_manual(values = colours, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(unique(df$stage)),
                       labels = unique(df$stage),
                       limits = c(1, length(unique(df$stage)))) +
    labs(x = "age", title = gene) +
    guides(fill = guide_legend(ncol = 2)) +
    ylab(glue("proportion {gene}+ cells")) +
    ylim(0, ymax) 
  
  return(p1)
  
}
