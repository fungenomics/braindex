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

# Vector specifying the order of clusters in the dendrogram
load("data/joint_mouse/ID_20190715_dendrogram_order.Rda")

# Load names of genes detected in mouse to provide choices in input
genes_mouse <- data.table::fread("data/joint_mouse/joint_mouse.gene_names.tsv", data.table = FALSE)$genes

# ---- Shiny settings ----

# Enable bookmarking
enableBookmarking(store = "url")