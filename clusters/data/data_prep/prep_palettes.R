
library(tidyr)
library(dplyr)
library(readr)

metadata <- data.table::fread("../joint_mouse/metadata_20190715.tsv", data.table = FALSE) %>% 
  select(Sample, Age, Species, Structure, Alias, Cell_type, Cluster, Colour, Cell_class, N_cells) %>% 
  mutate(Cluster_nounderscore = gsub("_", " ", Cluster)) %>% 
  write_tsv("../joint_mouse/metadata_20190715_select.tsv")

# Brain region palettes
pons_palette <- metadata %>%
  filter(Species == "Mouse" & grepl("[Pp]ons|[Hh]indbrain", Structure)) %>% 
  select(Cell_type, Colour) %>% 
  distinct() %>% 
  tibble::deframe()

cortex_palette <- metadata %>%
  filter(Species == "Mouse" & Structure == "Forebrain") %>%
  select(Cell_type, Colour) %>% 
  distinct() %>% 
  tibble::deframe()

x <- metadata %>% 
  select(Cluster, Colour)

joint_mouse_palette_refined <- x %>% 
  select(Cluster, Colour) %>% 
  tibble::deframe()

joint_mouse_palette_refined_nounderscore <- x %>% 
  select(Cluster, Colour) %>% 
  mutate(Cluster2 = gsub("_", " ", Cluster)) %>% 
  tibble::deframe()

rdbu <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100))

# Timepoint palettes
samples <- metadata %>% 
  filter(Species == "Mouse") %>% 
  distinct(Structure, Age) %>% 
  mutate(Sample = paste0(Structure, " ", Age)) %>% 
  pull(Sample)

ylgnbu <- RColorBrewer::brewer.pal(6, "YlGnBu") %>% tail(5)

timepoint_palette <- c(ylgnbu,
                       ylgnbu,
                       ylgnbu[1])

names(timepoint_palette) <- samples

save(pons_palette, cortex_palette, rdbu, timepoint_palette,
     joint_mouse_palette_refined,
     joint_mouse_palette_refined_nounderscore,
     file = "../joint_mouse/palettes.Rda")
