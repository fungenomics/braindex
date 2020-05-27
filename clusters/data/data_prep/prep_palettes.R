
library(tidyverse)

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

rdbu <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100))

save(pons_palette, cortex_palette, rdbu, file = "../joint_mouse/palettes.Rda")
