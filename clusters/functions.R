# Shiny helpers ----

# Custom functions ----

# TODO: write documentation for this function
get_embedding <- function(sample,
                          dr_cols,
                          cluster_column,
                          other_columns = NULL) {
  
  df <- feather::read_feather(glue("data/{sample}/{sample}.embedding_and_genes.feather"),
                              columns = c("Cell",
                                          dr_cols,
                                          cluster_column,
                                          other_columns))
  
  names(df)[4] <- "Cluster"
  
  return(df)
  
}

# TODO: write documentation for this function
get_expression <- function(sample,
                           embedding,
                           gene,
                           aggregate = TRUE) {
  
  # Add the gene expression levels to the embedding
  df <- cbind(embedding,
              feather::read_feather(glue("data/{sample}/{sample}.embedding_and_genes.feather"),
                                    columns = gene))
  
  if (aggregate) {
    
    if (length(gene) == 1) {
      
      df2 <- df[, 1:5]
      names(df2)[5] <- "Expression"
      
    } else if (length(gene) > 1) {
      
      # Take the mean of all the gene columns
      meanexp <- df[, 5:ncol(df)] %>% rowMeans
      
      df2 <- df
      df2$Expression <- meanexp
    
    }
    
    df <- df2
    
  }
  
  is_zero = FALSE
  
  non_zero_exp <- df %>% 
    filter(Expression != 0)
  
  # TRUE if non_zero_exp has no rows i.e. i.e. all zero expression
  if (dim(non_zero_exp)[1] == 0){ 
    
    is_zero = TRUE    
 
  }
  
  return(list(data = df, zero = is_zero))
  
}


#' Prepare input for bubble_plot
#' 
#' Load gene expression data from feather, tidy & optionally scale expression,
#' and return the dataframe required as input for the bubble_plot() function
#' 
#' @param gene Character vector, one or more genes of interest to plot
#' @param scale Logical, whether or not to linearly scale gene expression across
#' clusters to [0,1] to improve visualization. Default: TRUE
#' @param show_mean Logical, whether or not to display the mean expression of
#' given genes in a new bubble plot line. Default: FALSE
#' TODO: Use this to provide an option to download the underlying data.
#' 
#' @examples 
#' bubble_prep("Dlx1")
bubble_prep <- function(gene,
                        scale = TRUE,
                        show_mean = FALSE) {
  
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
  
  # Load the proportion of cells in each cluster in which each gene was detected,
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
    # However, letters take up more pixels than spaces, so do less padding
    # for genes with longer names
    mutate(Gene_padded = case_when(
      str_length(Gene) <= 5 ~ str_pad(Gene, 15, side = 'right', pad = " "),
      str_length(Gene) > 5 ~ str_pad(Gene, 12, side = 'right', pad = " ")
      )
    ) %>% 
    mutate(Gene_padded = factor(Gene_padded, levels = unique(.$Gene_padded))) %>% 
    
    # Order the clusters on the x-axis to match the dendrogram image
    mutate(Cluster = factor(Cluster, levels = dendrogram_order)) %>%
    
    filter(!is.na(Cluster)) %>% 
    
    # Convert NAs (undetected genes) to 0s -- this ensures all
    # clusters have a value for all genes, so that all clusters are plot,
    # even if the gene was undetected
    replace_na(list(Expression = 0, Pct1 = 0)) %>% 
    
    # Keep columns
    select(Gene, Cluster, Sample, Cell_type, Cell_class, N_cells, Expression, Pct1, Sample, Colour, Gene_padded)
  
  # Create & append set of rows containing mean expression over all selected genes
  if(show_mean) {
    
    # Create mean expression rows, preserving information for tooltip
    mean_exp <- df %>% 
      # We mean to group by cluster, but the other variables will be consistent for a given cluster
      group_by(Cluster, Sample, Cell_type, N_cells, Cell_class, Colour) %>%
      summarize(Expression = mean(Expression)) %>% 
      
      # Remove the Pct1 value from the mean expression
      # and label the mean expression
      mutate(Pct1 = 1, Gene = "MEAN", Gene_padded = "MEAN") 
    
    # Add the rows containing mean expression to the original dataframe,
    # removing duplicate rows and ordering them once more by user input,
    # except the mean which is placed at the bottom
    gene_order_padded <- levels(df$Gene_padded)
    df <- bind_rows(df, mean_exp) %>% 
      distinct(.) %>%
      mutate(Gene_padded = factor(Gene_padded, levels = c("MEAN", gene_order_padded)))
    
  }
  
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
#' @return A list containing a ggplot2 object and its legend (extracted with cowplot)
#'
#' @examples
#' bubble_prep("Dlx1") %>% bubbleplot()
#' 
#' @export
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
    scale_y_discrete(position = "right") #REMOVE GENE LABELS FROM PLOT - they will be separate
  
  gene_labels <- cowplot::get_y_axis(plot = p1, position = "right")
  
  p1 <- p1 + scale_y_discrete(labels = NULL)
  
  return(list(plot = p1, labels = gene_labels))
  
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
#' @param make_plotly Logical, whether or not to make an interactive (plotly)
#' version of the plot
#'
#' @return ggplot2 object or a plotly object, depending on make_plotly param
#'
#' @examples
#' ribbon_plot("Pdgfra", "joint_cortex")
ribbon_plot <- function(gene, 
                        region, 
                        ymax = NA, 
                        make_plotly = FALSE) {
  
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
    distinct(Age, Cell_type, frac, total) %>% 
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
  
  df <- left_join(df, select(ribbon_df_celltype_frac, cluster = Cell_type, stage = Age, frac, total)) %>% 
    # Complete cases when genes were not detected in certain timepoints/clusters
    # by replacing with a zero
    mutate(frac = replace_na(frac, 0)) %>% 
    left_join(select(ribbon_df_cum_frac, stage = Age, cumfrac))
  
  df$xpos = match(df$stage, unique(timepoints2))
  
  p1 <- df %>%
    # Need to specify group or the text attribute with glue causes errors
    ggplot(aes(x = xpos, y = frac, fill = cluster, group = cluster,
               text = glue("{total*frac} {gene}+ cells out of {total} cells at this time point"))) +
    geom_area(stat = "identity") +
    scale_fill_manual(values = colours, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(unique(df$stage)),
                       labels = unique(df$stage),
                       limits = c(1, length(unique(df$stage)))) +
    labs(x = "age", title = gene) +
    guides(fill = guide_legend(ncol = 2)) +
    ylab(glue("proportion {gene}+ cells")) +
    ylim(0, ymax) 
  
  if(make_plotly) {
    return ((ggplotly(p1,
                     # Display cluster (group) and info on number of cells (text) in tooltips
                     tooltip = c("group", "text")) %>%
              
              # Add hovers both on points as well as filled areas of the plot
              # Changing it to hoveron="fills" only causes a known issue, see:
              # https://github.com/ropensci/plotly/issues/1641 
              style(hoveron="points+fills") 
            ) %>% 
              
              # Customize the modebar on the plotly object to hide certain buttons, 
              # remove the plotly logo, and toggle spike lines on by default
              config(modeBarButtonsToRemove = c("hoverCompareCartesian", 
                                                "hoverClosestCartesian",
                                                "toImage"),
                     displaylogo = FALSE) %>% 
              layout(yaxis = list(showspikes = FALSE,
                                  spikethickness = 1.5,
                                  spikedash = "solid"),
                     xaxis = list(showspikes = FALSE,
                                  spikethickness = 1.5,
                                  spikedash = "solid")))
  }  else {
    return(p1)
  }
  
}

#TODO: finish documentation for this function
#' Generate dimensionality reduction plot from data embedding
#' 
#' @param embedding ...
#' @param colour_by String, variable to colour the plot by. Default: NULL
#' @param colours Character vector, colour palette to use for plot. Default: NULL
#' @param colour_by_type String, colour palette type, either "continuous" or 
#' "discrete". Default: "discrete"
#' @param label Logical, whether or not to label clusters in plot. Default: TRUE
#' @param point_size Numeric, size of points in mm. Default: 0.4
#' @param alpha Numeric, transparency of points (from 0 to 1). Default: 0.8
#' @param legend Logical, whether or not to include a legend in the plot.
#' Default: FALSE if label = TRUE and colour_by = NULL
#' NOT USED? @param label_repel Logical, ... Default: TRUE
#' @param label_size Numeric, font size of cluster labels. Default: 4
#' @param cells Default: NULL
#' @param order_by Default: Null
#' @param clusters_to_label Default: NULL
#' @param hide_ticks Logical, whether or not to hide axis ticks. Default: TRUE
#' @param title String, title to be displayed on the plot. Default: NULL
#' NOTE USED? @param label_short Logical, ... Default: FALSE
#' @param na_color String, colour of NA values in plot. Default: "gray80"
#' @param limits Numeric vector, limits used for a continuous colour palette. 
#' Default: NULL
#' @param hide_axes Logical, whether or not to hide the plot axes. Default: FALSE
#' @param show_n_cells Logical, ... Default: FALSE
#' 
#' @return A list containing a ggplot object and a list of cluster centers
#' 
#' @export
dr_plot <- function(embedding,
                    colour_by = NULL,
                    colours = NULL,
                    colour_by_type = "discrete",
                    label = TRUE,
                    point_size = 0.4,
                    alpha = 0.8,
                    legend = ifelse((is.null(colour_by)) && (label), FALSE, TRUE),
                    label_repel = TRUE,
                    label_size = 4,
                    cells = NULL,
                    order_by = NULL,
                    clusters_to_label = NULL,
                    hide_ticks = TRUE,
                    title = NULL,
                    label_short = FALSE,
                    na_color = "gray80",
                    limits = NULL,
                    hide_axes = FALSE,
                    show_n_cells = FALSE) {
  
  order_by <- colour_by
  
  dr_cols <- colnames(embedding)[c(2, 3)] %>% 
    # Remove underscores
  {gsub("_", " ", .)}
  
  colnames(embedding)[c(2, 3)] <- c("dim1", "dim2")
  
  if (!is.null(cells)) embedding <- embedding %>% filter(Cell %in% cells)
  if (!is.null(order_by)) embedding <- embedding %>% arrange_(order_by)
  
  if (show_n_cells) {
    
    embedding <- embedding %>%
      group_by(Cluster) %>%
      mutate(n_cells = n()) %>% 
      mutate(Cluster2 = paste0(Cluster, " (", n_cells, ")"))
    
    # Rename the palette
    names(colours) <- paste0(names(colours), " (", 
                             embedding %>% distinct(cluster, n_cells) %>% arrange(desc(n_cells)) %>% .$n_cells,
                             ")")
    
  } else {
    
    embedding$Cluster2 <- embedding$Cluster
    
  }
  
  gg <- ggplot(embedding, aes(x = dim1, y = dim2))
  
  # Deal with the palette
  if (is.null(colour_by)) {
    
    gg <- gg +
      geom_point(aes(colour = Cluster2), size = point_size, alpha = alpha) +
      scale_color_manual(values = colours)
    
  } else {
    
    if (is.null(limits)) lims <- c(NA, NA) # NAs refer to the current min & max values
    else lims <- limits
    
    gg <- gg +
      geom_point(aes_string(colour = colour_by), size = point_size, alpha = alpha)
    
    if (!is.null(colours)) {
      
      if (colour_by_type == "discrete") gg <- gg + scale_color_manual(values = colours, na.value = na_color)
      else if (colour_by_type == "continuous") {
        
        gg <- gg + scale_color_gradientn(colours = colours,
                                         na.value = na_color,
                                         limits = lims)
      }
      
    } else {
      
      if (colour_by_type == "continuous") { # Otherwise for discrete, default ggplot2 colours are used
        
        gg <- gg + scale_color_gradientn(colours = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                         na.value = na_color,
                                         limits = lims)
        
      }
    }
  }
  
  # Store the center points (medians) of each cluster
  centers <- embedding %>%
    group_by(Cluster) %>%
    summarise(center_x = median(dim1),
              center_y = median(dim2))
  
  # Label clusters at cluster centers
  if (label) {
    
    gg <- gg + ggrepel::geom_label_repel(data = centers,
                                         aes(x = center_x, y = center_y),
                                         label = centers$Cluster,
                                         size = label_size,
                                         segment.color = 'grey50',
                                         fontface = 'bold',
                                         alpha = 0.8,
                                         segment.alpha = 0.8,
                                         label.size = NA,
                                         force = 2,
                                         segment.size = 0.5,
                                         arrow = arrow(length = unit(0.01, 'npc')))
    
  }
  
  gg <- gg + theme_min()
  
  # Set the right axes titles
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  else gg <- gg + 
    labs(x = dr_cols[1],
         y = dr_cols[2],
         colour = "Cluster")
  
  if (hide_ticks) gg <- gg + noTicks()
  
  # More aesthetics
  if (!legend) gg <- gg + theme(legend.position = "none")
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  
  # Return the cluster centers for hover functionality on plot
  return(list("plot" = gg, "centers" = centers))
  
}

#' Colour cells in t-SNE or PCA space by gene expression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to t-SNE space, but see the \code{reduction} argument for
#' how to plot in PCA space instead. This function is based on \code{Seurat::FeaturePlot}.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param genes String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list seurat@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#' @param dim1 Numeric, index of dimension from \code{reduction} to plot on
#' the x-axis. e.g. to plot the 3rd PC on the x-axis, pass 3. Default: 1.
#' @param dim2 Numeric, like \code{dim2}, but for the y-axis. Default: 2.
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scale. Default: redgrey.
#' @param title (Optional) String specifying the plot title
#' @param alpha Numeric, fixed alpha for points. Default: 0.6
#' @param point_size Numeric, size of points in scatterplot. Default: 1. (A smaller
#' value around 0.5 is better for plots which will be viewed at small scale.)
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param legend Logical, whether or not to plot legend. Default: TRUE
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: TRUE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param limits (Optional) A numeric vector of length two providing the limits to
#' use for the colour scale (documentation
#' from \code{\link[ggplot2]{continous_scale}}. Default: 0 and max of the data.
#'
#' @export
#' @return A ggplot object
feature_plot <- function(df,
                         label = TRUE,
                         palette = "redgrey",
                         title = NULL,
                         alpha = 0.6,
                         label_repel = TRUE,
                         label_size = 4,
                         legend = TRUE,
                         hide_ticks = FALSE,
                         hide_axes = FALSE,
                         y_max = NULL,
                         label_short = FALSE,
                         return_df = FALSE,
                         point_size = 1,
                         order_size = NULL,
                         legend_title = "Expression") {
  
  dr_cols <- colnames(df)[c(2, 3)] %>% 
    # Remove underscores
  {gsub("_", " ", .)}
  
  colnames(df)[c(2, 3)] <- c("dim1", "dim2")
  
  df <- df %>% dplyr::arrange(Expression)
  
  # Set limits: if not provided, use default min/max
  if (is.null(y_max)) limits <- c(NA, NA)
  else {
    
    limits <- c(0, y_max)
    
  }
  
  # If ALL gene expression values = 0, change all expression values to NA
  # to control their colour via the na.value parameter of ggplot
  non_zero_exp <- df %>% 
    filter(Expression != 0)
  
  if (dim(non_zero_exp)[1] == 0 && # TRUE if non_zero_exp has no rows
      palette != "viridis"){ # Don't set NA for viridis palette
    
    df <- df %>% 
      mutate(Expression = as.numeric(NA))
    
  }
  
  # Plot using the palette chosen by the user
  gg <- df %>%
    ggplot(aes(x = dim1, y = dim2)) +
    geom_point(aes(colour = Expression), size = point_size, alpha = alpha)
  
  if (palette == "viridis") {
    
    gg <- gg + viridis::scale_color_viridis(limits = limits)
    
  } else if (palette == "blues") {
    
    gg <- gg + scale_colour_gradientn(
      colours = RColorBrewer::brewer.pal(n = 8, name = "Blues"),
      limits = limits,
      # NA values set to the colour assigned to 0 in this palette
      na.value = "#F7FBFF") 
    
  } else if (palette == "redgrey") {
    
    # NOTE: palette chosen is not the default gradient from gray -> red
    # but sets a midpoint at a lighter colour
    gg <- gg + scale_color_gradientn(
      colours = grDevices::colorRampPalette(colors = c("gray83", "#E09797", "red"))(n = 200),
      limits = limits,
      # NA values set to the colour assigned to 0 in this palette
      na.value = "grey83") 
    
  } else if (palette == "rdbu") {
    
    gg <- gg + scale_color_gradientn(
      colours = tail(rdbu, 70),
      limits = limits,
      # NA values set to the colour assigned to 0 in this palette
      na.value = "#99C8E0") 
    
  }
  
  # If all expression values equaled 0, they were set to NA, and 
  # the colour palette legend disappeared.
  # Add caption below the plot to clarify that all expression = 0 
  if (dim(non_zero_exp)[1] == 0){ # TRUE if non_zero_exp has no rows

    gg <- gg + labs(
      # Newlines (\n) to fill same space as legend, for plot alignment purposes
      caption = "Gene expression = 0 for all cells in this plot. \n \n \n \n")
    
  }
  
  # Label clusters
  if (label) {
    
    centers <- df %>%
      group_by(cluster) %>%
      summarise(mean_x = median(dim1),
                mean_y = median(dim2))
    gg <- gg + ggrepel::geom_label_repel(data = centers,
                                         aes(x = mean_x, y = mean_y),
                                         label = centers$cluster,
                                         size = label_size,
                                         segment.color = 'grey50',
                                         fontface = 'bold',
                                         alpha = 0.8,
                                         segment.alpha = 0.8,
                                         label.size = NA,
                                         force = 2,
                                         segment.size = 0.5,
                                         arrow = arrow(length = unit(0.01, 'npc')))
    
  }
  
  gg <- gg + theme_min()
  
  # Set the right axes titles
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  else gg <- gg + xlab(dr_cols[1]) + ylab(dr_cols[2])
  
  if (hide_ticks) gg <- gg + noTicks()
  
  # More aesthetics
  if (!legend) gg <- gg + theme(legend.position = "none")
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  
  gg <- gg + labs(colour = legend_title)
  
  return(gg)
  
}

#' Generate a violin plot of single cell data
#' 
#' @param df Dataframe, contains data shown in plot
#' @param palette Character vector, colour palette to be used for the plot
#' @param scale String, determining the scale input for geom_violin, i.e. 
#' the parameter that will remain the same between violins. Default: "width",
#' other possible values are "area" and "count"
#' @param points Logical, whether or not to show points in plot. Default: FALSE
#' @param point_size Numeric, indicating the size of points in mm. Default: 0.4
#' @param y_lab String, label for the y-axis of the plot. Default: "Normalized 
#' expression"
#' 
#' @return a ggplot object
vln <- function(df,
                palette,
                scale = "width",
                points = FALSE,
                point_size = 0.4,
                y_lab = "Normalized expression") {
  
  cluster_order <- df %>% 
    group_by(Cluster) %>% 
    summarize(Mean_exp = mean(Expression)) %>% 
    arrange(desc(Mean_exp)) %>% 
    pull(Cluster)
  
  gg <- df %>% 
    filter(!grepl("BLACKLIST", Cluster)) %>% 
    mutate(Cluster = factor(Cluster, levels = cluster_order)) %>% 
    ggplot(aes(x = Cluster, y = Expression)) +
    geom_violin(scale = scale, aes(fill = Cluster)) +
    labs(x = "Cluster",
         y = y_lab,
         fill = "Cluster") +
    scale_fill_manual(values = palette) +
    theme_min() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  
  if (points) gg <- gg + geom_jitter(size = point_size)
  
  return(gg)
  
}

#' Remove ticks from the axis of a ggplot
#' 
#' @example 
#' gg <- gg + noTicks() 
noTicks <- function() {
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
}

#' Determine if a background colour is dark enough to warrant white text
#' 
#' @param hex_color String, colour in hex colour format e.g. #000000
#' 
#' @return TRUE if the colour is dark enough (arbitrary)
dark <- function(hex_color) {
  
  red <- substr(hex_color, 2, 2)
  green <- substr(hex_color, 4, 4)
  blue <- substr(hex_color, 6, 6)
  dark_nums <- c(0:8)
  
  if ((red %in% dark_nums && blue %in% dark_nums) || 
      (red %in% dark_nums && green %in% dark_nums) ||
      (green %in% dark_nums && blue %in% dark_nums)) {
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
    
  }
}

#' Add ticks below a bar plot to categorize x axis into less granular categories
#' 
#' @param df Dataframe, containing the data to use
#' [...]
#' 
#' @example 
#' plot + add_class_ticks(df, unique(df$Cell_class), palette = palettes$Cell_class,
#'                        start = -50, sep = 100, height = 500, label_x_pos = -9, fontsize = 3.5)
#'
add_class_ticks <- function(df, classes, height, sep, start, label_x_pos, palette = NULL, fontsize = 3) {
  
  # Set up our limits
  n <- length(classes)
  tops     <- seq(start, by = - (height + sep), length.out = n)
  bottoms  <- seq(start - height, by = - (height + sep), length.out = n)
  mids     <- map2_dbl(tops, bottoms, ~ mean(c(.x, .y)))
  betweens <- seq(start - (height + sep/2), by = - (height + sep), length.out = n - 1)
  
  if (is.null(palette)) palette <- rep("black", n)
  
  # Make a dataframe for tick positions
  df$y_top <- NA
  df$y_bottom <- NA
  
  for (i in seq_along(classes)) {
    
    df[df$Cell_class == classes[i], ]$y_top <- tops[i]
    df[df$Cell_class == classes[i], ]$y_bottom <- bottoms[i]
    
  }
  
  # Make a dataframe for class labels
  df2 <- data.frame(Class = classes,
                    x = label_x_pos,
                    y = mids)
  
  # Adding ggplot2 elements together
  # https://stackoverflow.com/questions/56405904/how-to-add-ggproto-objects-together-and-save-for-later-without-call-to-ggplot
  list(geom_segment(data = df,
                    mapping = aes(x = Cluster, y = y_top,
                                  xend = Cluster, yend = y_bottom),
                    size = 1,
                    colour = "gray50"),
       geom_hline(yintercept = 0, colour = "gray90"),
       geom_hline(yintercept = betweens, linetype = "dotted", size = 0.4, colour = "gray60"),
       geom_text(data = df2, mapping = aes(x = x, y = y, label = Class, colour = Class), size = fontsize, fontface = "bold", hjust = "left"),
       scale_colour_manual(values = palette))
  
}
