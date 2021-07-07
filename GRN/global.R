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
#library(rcytoscapejs2) # downloaded from https://github.com/uc-bd2k/rcytoscapejs2
#could probably get rid of rcytoscape
library(glue)
library(GGally)
library(sna)
library(igraph)
library(intergraph)
library(purrr)
library(shinyWidgets)

# Data
load("data/joint_cortex/cortex_prep.Rda") # a list, data_cortex
load("data/joint_pons/pons_prep.Rda")     # a list, data_pons
load("data/shared/common_prep.Rda") # metadata and colour_palettes

data_ct_e12 <- readRDS("data/ct_e12/ct_e12_prep.Rds") # a list, data_ct_e12
data_ct_e15 <- readRDS("data/ct_e15/ct_e15_prep.Rds")
data_ct_p0 <- readRDS("data/ct_p0/ct_p0_prep.Rds")
data_ct_p3 <- readRDS("data/ct_p3/ct_p3_prep.Rds")
data_ct_p6 <- readRDS("data/ct_p6/ct_p6_prep.Rds")

data_po_e12 <- readRDS("data/po_e12/po_e12_prep.Rds") 
data_po_e15 <- readRDS("data/po_e15/po_e15_prep.Rds")
data_po_p0 <- readRDS("data/po_p0/po_p0_prep.Rds")
data_po_p3 <- readRDS("data/po_p3/po_p3_prep.Rds")
data_po_p6 <- readRDS("data/po_p6/po_p6_prep.Rds")

# Custom functions
source("functions.R")

#allows server side save states
enableBookmarking(store = "url")
