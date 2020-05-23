
load("data/shiny_precomputed_resources.Rda")

plot_gene_in_time <- function(gene,
                              span = 1,
                              points = TRUE,
                              trend = TRUE,
                              min_prop = 0.1,
                              x_axis = "even",
                              jitter_width = 0.3,
                              facet_cell_type = FALSE,
                              point_size = 1) {
  
	set.seed(100)

  y <- feather::read_feather("data/joint_cortex2.feather", c("orig.ident", "Cell", gene))
  
  if (length(gene) > 1) {
    
    mean_expr <- rowMeans(y[, 3:ncol(y)])
    y$Gene <- mean_expr
    
  } else  colnames(y)[3] <- "Gene"
  
  y$Cell_type <- jc_all_meta_labelled$Cell_type
  y$Colour <- jc_all_meta_labelled$Colour
  y$Old_cluster <- jc_all_meta_labelled$ID_20190715_with_blacklist
  
  df <- y %>%
    # Filter clusters
    filter(!grepl("BLACKLIST", Old_cluster)) %>% 
    filter(!grepl("Mixed|Ery|Unresolved|Split|Active|Cck", Cell_type)) %>% 
    mutate(orig.ident2 = orig.ident) %>%
    separate(orig.ident2, into = c("region", "age"), sep = " ") %>% 
    left_join(timepoints, by = "age") %>% 
    # Calculate proportions
    group_by(Cell_type) %>% 
    mutate(n = n()) %>% 
    mutate(Prop_detected = sum(Gene >  0)/n) %>% 
    ungroup() %>% 
    arrange(desc(Prop_detected)) %>% 
    mutate(Group = paste0(Cell_type, " (", round(Prop_detected, 2)*100, "%)")) %>% 
    mutate(Group = factor(Group, levels = unique(.$Group))) %>% 
    filter(Prop_detected >= min_prop)
  
  if (x_axis == "even") {
    
    df$X <- df$numeric
    breaks <- timepoints$numeric
    xlim <- c(1, 9)
    
  }  else if (x_axis == "representative") {
    
    df$X <- df$numeric3
    breaks <- timepoints$numeric3
    xlim <- c(-11, 6)
    
  }
  
  palette_with_prop <- df %>%
    select(Group, Colour) %>% 
    distinct(Group, .keep_all = TRUE) %>% 
    tibble::deframe()
  
  p1 <- df %>%
    ggplot(aes(x = X, y = Gene))
  
  if (points) p1 <- p1 +
    geom_jitter(aes(colour = Group),
                alpha = 0.5,
                width = jitter_width,
                size = point_size)
  
  if (trend) p1 <- p1 + 
    stat_smooth(geom ='line', data = df %>% filter(Gene > 0),
                mapping = aes(x = X, y = Gene, group = Group, alpha = Prop_detected, colour = Group),
                se = FALSE,
                method = "loess",
                size = 2,
                span = span)
  
  p1 <- p1 +
    scale_colour_manual(values = palette_with_prop, name = "Cell type") + 
    scale_alpha_continuous(range = c(0.4, 1)) +
    guides(alpha = FALSE) +
    ylab(ifelse(length(gene) == 1, "Expression", "Mean expression")) +
    xlab("Age") +
    ggtitle(ifelse(length(gene) == 1, gene, "Mean expression")) +
    ylim(c(0, 4)) +
    scale_x_continuous(breaks = breaks, labels = timepoints$age, limits = xlim) +
    theme(axis.title   = element_text(size = 12, face = "bold"),
          axis.text    = element_text(size = 10),
          title        = element_text(size = 12, face = "bold.italic"),
          legend.title = element_text(size = 12, face = "bold"))
  
  if (facet_cell_type) p1 <- p1 + facet_wrap(~ Group) + theme(legend.position = "none")
  
  legend <- cowplot::get_legend(p1)
  
  return(list("plot" = p1 + theme(legend.position = "none"),
              "legend" = legend))
  
  
}



palette_group <- c("Excitatory neuron"  = "#c9110e",
                   "Inhibitory neuron"  = "#4a1cc9",
                   "Oligodendrocyte"    = "#b7dd5f",
                   "Astro-ependymal"    = "#0b7c68")

plot_gene_in_pseudotime <- function(gene,
                                    span = 0.7,
                                    points = TRUE,
                                    trend = TRUE,
                                    min_prop = 0.1,
                                    facet_cell_type = FALSE,
                                    point_size = 1) {
 
  set.seed(100)

  y <- feather::read_feather("data/joint_cortex.embedding_and_genes.feather",
                             c("orig.ident", "Cell", gene))
  
  if (length(gene) > 1) {
    
    mean_expr <- rowMeans(y[, 3:ncol(y)])
    y$Gene <- mean_expr
    
  } else  colnames(y)[3] <- "Gene"
  
  df <- bind_cols(y, jc_meta_pseudotime) %>% 
    filter(!is.na(pseudotime))
  
  df <- df %>%
    filter(!grepl("BLACKLIST", Cell_type)) %>% 
    filter(!grepl("Mixed|Ery|Unresolved|Split|Active|Cck|hem|Cajal|Menin|Endo|Stria|Thal|Peri|neuronal", Cell_type)) %>% 
    # Calculate proportions
    group_by(Lineage) %>% 
    mutate(n = n()) %>% 
    mutate(Prop_detected = sum(Gene >  0)/n) %>% 
    ungroup() %>% 
    arrange(desc(Prop_detected)) %>% 
    mutate(Lineage = factor(Lineage, levels = unique(.$Lineage))) %>% 
    filter(Prop_detected >= min_prop)
  
  p1 <- df %>%
    ggplot(aes(x = pseudotime_scaled, y = Gene))
  
  if (points) p1 <- p1 +
    geom_point(aes(colour = Cell_type),
               alpha = 0.5,
               size = point_size)
  
  if (trend) p1 <- p1 + 
    stat_smooth(geom ='line', data = df,
                mapping = aes(x = pseudotime_scaled, y = Gene, group = Lineage, alpha = Prop_detected, colour = Lineage),
                se = FALSE,
                method = "loess",
                size = 2,
                span = span)
  
  p1 <- p1 +
    scale_colour_manual(values = c(palette_group, palette_cortex), name = "Cell type") + 
    scale_alpha_continuous(range = c(0.5, 1)) +
    guides(alpha = FALSE) +
    ylab(ifelse(length(gene) == 1, "Expression", "Mean expression")) +
    xlab("Pseudotime") +
    ggtitle(ifelse(length(gene) == 1, gene, "Mean expression")) +
    ylim(c(0, 4)) +
    scale_x_continuous(breaks = c(10, 90), labels = c("Early", "Late")) +
    theme(axis.title   = element_text(size = 12, face = "bold"),
          axis.text    = element_text(size = 10),
          title        = element_text(size = 12, face = "bold.italic"),
          legend.title = element_text(size = 12, face = "bold"))
  
  if (facet_cell_type) p1 <- p1 + facet_wrap(~ Lineage) + theme(legend.position = "none")
 
  legend <- cowplot::get_legend(p1)

  return(list("plot" = p1 + theme(legend.position = "none"),
              "legend" = legend))
  
}



ribbon_plot <- function(gene, mode = 1, ymax = NA) {
  
  breaks <- timepoints$numeric
  xlim <- c(1, 9)

  ribbon_df <- feather::read_feather("data/joint_cortex2.feather",
                                     c("orig.ident", "Cell", gene)) %>% 
    rowwise() %>% 
    mutate(Cell = stringr::str_extract(Cell, "[ACGT]{16}")) %>% 
    ungroup()
  
  ribbon_df$Cell_type <- jc_all_meta_labelled$Cell_type
  ribbon_df$Colour <- jc_all_meta_labelled$Colour
  
  ribbon_df <- ribbon_df %>% 
    filter(!grepl("Mixed|Ery|Unresolved|Split|Active|Cck", Cell_type)) %>% 
    mutate(orig.ident = factor(orig.ident, levels = c("Forebrain E10.5",
                                                      "Forebrain E12.5",
                                                      "Forebrain E13.5",
                                                      "Forebrain E15.5",
                                                      "Forebrain E16.5",
                                                      "Forebrain E18.5",
                                                      "Forebrain P0",
                                                      "Forebrain P3",
                                                      "Forebrain P6"))) %>% 
    arrange(orig.ident)
  
  ribbon_df$gene <- ribbon_df[[gene]]
  
  if (mode == 1) {
    
    ribbon_df_celltype_frac <- ribbon_df %>% 
      group_by(orig.ident) %>% 
      mutate(total = n()) %>% 
      group_by(orig.ident, Cell_type) %>% 
      mutate(frac = sum(gene > 0) / total) %>% 
      distinct(orig.ident, Cell_type, frac) %>% 
      ungroup()
    
    ribbon_df_cum_frac <- ribbon_df %>% 
      group_by(orig.ident) %>% 
      summarize(cumfrac = sum(gene > 0) / n()) %>% 
      ungroup()
    
  } else if (mode == 2) {
    
    ribbon_df_celltype_frac <- ribbon_df %>% 
      filter(gene > 0) %>% 
      group_by(orig.ident) %>% 
      mutate(total = n()) %>% 
      group_by(orig.ident, Cell_type) %>% 
      mutate(frac = n() / total) %>% 
      distinct(orig.ident, Cell_type, frac) %>% 
      ungroup()
    
    ribbon_df_cum_frac <- data.frame(orig.ident = unique(ribbon_df$orig.ident),
                                     cumfrac = 1)
    
  }
  
  timepoints2 <- ribbon_df$orig.ident
  clusters <- ribbon_df$Cell_type
  colours <- palette_cortex
  
  df = data.frame(cluster = rep(unique(clusters), length(unique(timepoints2))),
                  stage = do.call(c, lapply(as.character(unique(timepoints2)), rep, times = length(unique(clusters)))))
  
  df$ranking = match(df$cluster, names(colours))
  df = df[order(df$stage, df$ranking),]
  
  df <- left_join(df, select(ribbon_df_celltype_frac, cluster = Cell_type, stage = orig.ident, frac)) %>% 
    mutate(frac = replace_na(frac, 0)) %>% 
    left_join(select(ribbon_df_cum_frac, stage = orig.ident, cumfrac))
  
  df$xpos = match(df$stage, unique(timepoints2))
  
  p1 <- df %>%
    ggplot(aes(x = xpos, y = frac, fill = cluster)) +
    geom_area(stat = "identity") +
    scale_fill_manual(values = colours, drop = FALSE, name = "") +
    scale_x_continuous(breaks = breaks, labels = timepoints$age, limits = xlim) +
    labs(x = "age", title = gene) +
    guides(fill = guide_legend(ncol = 3))
  
  
  if (mode == 1) {
    
    p1 <- p1  + ylab(glue("proportion {gene}+ cells")) +
      ylim(0, ymax)
    
  } else if (mode == 2) {
    
    p1 <- p1 + ylab(glue("breakdown of {gene}+ cells"))
    
  }
  
  legend <- cowplot::get_legend(p1)
  
  return(list("plot" = p1 + theme(legend.position = "none"),
              "legend" = legend))

  
}
