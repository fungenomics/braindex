# This file will be run first, and anything in this file will be available
# for the other R files required for the shiny app (e.g. server.R, ui.R)

# ---- Load required packages ----
library(feather)
library(tidyr)
library(dplyr)
library(cowplot)
library(glue)
library(stringr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(readr)
library(shinyWidgets)
library(plotly)
library(shinycssloaders)
library(shiny)
library(reactable)

# ---- Set-up / load common data ----

# Cluster-level metadata
metadata <- data.table::fread("data/joint_mouse/metadata_20190715_select.tsv",
                              data.table = FALSE)

# Red-blue gradient and brain region colour palettes
load("data/joint_mouse/palettes.Rda")

cortex_palette_joint <- readRDS("data/joint_cortex/joint_cortex.palette_ID_20190715_joint_clustering.Rds")
pons_palette_joint   <- readRDS("data/joint_pons/joint_pons.palette_ID_20190715_joint_clustering.Rds")

# Joint mouse colour palette
load("data/joint_mouse/joint_mouse.palette_ID_20190715.Rda")

# General cell type palette
general_palette <- c("Progenitors/cyc." = "#ffaf49",
                       "Oligodendrocytes" = "#b7dd5f",
                       "Astrocytes" = "#00a385",
                       "Ependymal" = "#8ee5cf",
                       "Neurons" = "#840200",
                       "Non-neuroect." = "gray40",
                       "Other" = "gray60")

# Vector specifying the order of clusters in the dendrogram
load("data/joint_mouse/ID_20190715_dendrogram_order.Rda")

# Load names of genes detected in mouse - genes for which there is data in atlas
genes_mouse <- data.table::fread("data/joint_mouse/joint_mouse.gene_names.tsv", data.table = FALSE)$genes

# Load all genes in mouse annotation - to validate input from users & provide as choices
# Some of these genes may not have corresponding data in the atlas - 
# i.e. genes_mouse (above) is a subset of genes_anno
genes_anno <- data.table::fread("data/all_mm10_genes.txt", header = FALSE, data.table=FALSE)
names(genes_anno) <- "Genes"
genes_anno <- genes_anno[['Genes']]

# ---- Shiny settings ----

# -- Enable bookmarking
enableBookmarking(store = "url")

# ---- R settings ----

# -- Overwrite pheatmap function so the heatmap labels can flipped and 45 degrees
library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with a new version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))