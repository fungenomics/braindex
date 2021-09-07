# Libraries
library(rintrojs)
library(shiny)
library(shinyjs)
library(tidyr)
library(stringr)
library(tibble)
library(dplyr)
library(readr)
library(feather)
library(plotly)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(DT)
library(glue)
library(GGally)
library(sna)
library(igraph)
library(intergraph)
library(purrr)
library(shinyWidgets)

# Data
load("data/joint_cortex_extended/cortex_extended_prep.Rda") # a list, data_cortex_extended
load("data/joint_pons_extended/pons_extended_prep.Rda")     # a list, data_pons_extended
load("data/shared/common_prep.Rda") # metadata and colour_palettes

#a list of R objects for each per-sample SCENIC dataset prepared in data_prep.R
#Cortex
data_ct_e10 <- readRDS("data/ct_e10/ct_e10_prep.Rds")
data_ct_e12 <- readRDS("data/ct_e12/ct_e12_prep.Rds") 
data_ct_e13 <- readRDS("data/ct_e13/ct_e13_prep.Rds")
data_ct_e15 <- readRDS("data/ct_e15/ct_e15_prep.Rds")
data_ct_e16 <- readRDS("data/ct_e16/ct_e16_prep.Rds")
data_ct_e18 <- readRDS("data/ct_e18/ct_e18_prep.Rds")
data_ct_p0 <- readRDS("data/ct_p0/ct_p0_prep.Rds")
data_ct_p3 <- readRDS("data/ct_p3/ct_p3_prep.Rds")
data_ct_p6 <- readRDS("data/ct_p6/ct_p6_prep.Rds")

#Pons
data_po_e10 <- readRDS("data/po_e10/po_e10_prep.Rds")
data_po_e12 <- readRDS("data/po_e12/po_e12_prep.Rds")
data_po_e13 <- readRDS("data/po_e13/po_e13_prep.Rds")
data_po_e15 <- readRDS("data/po_e15/po_e15_prep.Rds")
data_po_e16 <- readRDS("data/po_e16/po_e16_prep.Rds")
data_po_e18 <- readRDS("data/po_e18/po_e18_prep.Rds")
data_po_p0 <- readRDS("data/po_p0/po_p0_prep.Rds")
data_po_p3 <- readRDS("data/po_p3/po_p3_prep.Rds")
data_po_p6 <- readRDS("data/po_p6/po_p6_prep.Rds")


#order of clusters in the dendrograms 
dend_order_joint_cortex_extended <- readRDS("data/joint_cortex_extended/dendrogram_order_joint_extended_forebrain.Rds")
dend_order_joint_pons_extended <- readRDS("data/joint_pons_extended/dendrogram_order_joint_extended_pons.Rds")

#load ribbon_plots of proportion of cells across developmental time in pons and forebrain
load("data/shared/timeseries_proportion_plots.Rda")

# Custom functions
source("functions.R")

#allows URL save states
enableBookmarking(store = "url")

#this vector is used repeatedly in app as the selection options for per timepoint data
#visualisation, assigned here to make it accessible to all of app
dev_time_points <- c("e10", "e12", "e13", "e15", "e16", "e18", "p0", "p3", "p6")

#overwrite pheatmap function that draws the column names on the heatmap
#so that the heatmap labels can be 45 degrees
library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

#palette for bubble plot
rdbu <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100))
