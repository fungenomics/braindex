# Libraries
library(rintrojs)
library(shiny)
#library(shinyjs)
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
library(rcytoscapejs2) # downloaded from https://github.com/uc-bd2k/rcytoscapejs2
#could probably get rid of rcytoscape
library(glue)
library(GGally)
library(sna)
library(igraph)
library(integraph)
library(purrr)

# Data
load("data/joint_cortex/cortex_prep.Rda") # a list, data_cortex
load("data/joint_pons/pons_prep.Rda")     # a list, data_pons
load("data/shared/common_prep.Rda") # metadata and colour_palettes

# Custom functions
source("functions.R")