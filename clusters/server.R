
# For the plot cache
# NOTE:
# To get this to work, need to:
# 1) Make the cache folder, i.e. mkdir cache
# 2) On the server, give the shiny user permissions, i.e. chmod -R a=rwx cache
shinyOptions(cache = diskCache("./cache"))

library(cowplot)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("functions.R")
source("../style.R")

# Set a default ggplot2 theme for the app, from cowplot
ggplot2::theme_set(theme_min())

server <- function(input, output, session) {
  
  
  
  # Capture all input from this tab as a list in case we want to add
  # more options in the future
  input_new <- eventReactive(input$update, {
    
    # Inputs to use as is
    l <- list(
      "gene"   = input$gene,
      "scale"  = input$bubble_scale,
      "size"   = input$bubble_size,
      "region" = input$region,
      "label_clusters" = input$label_clusters,
      "ft_palette"     = input$feature_palette
    )
    
    # Get the columns for the appropriate type of dim red
    if      (input$dr == "tSNE") l$dr <- c("tSNE_1", "tSNE_2")
    else if (input$dr == "UMAP") l$dr <- c("UMAP1", "UMAP2")
    else if (input$dr == "PCA")  l$dr <- c("PC1", "PC2")
    
    # Get the clustering to use for the joint analysis tab
    # Option 1: Clustering done at the joint analysis level
    # Option 2: Clustering done per sample (for the forebrain, there was some refinement done so the latest
    # version of labels is different than in the pons)
    if (input$dr_clustering == "Clustering at the region level") {
      
      l$clust <- "ID_20190715_joint_clustering"
      
      if       (input$region == "joint_cortex") l$clust_palette <- cortex_palette_joint
      else if (input$region == "joint_pons")    l$clust_palette <- pons_palette_joint
      
    } else {
      
      l$clust_palette <- joint_mouse_palette_refined
      
      if      (input$region == "joint_cortex") l$clust <- "ID_20190730_with_blacklist_and_refined"
      else if (input$region == "joint_pons")   l$clust <- "ID_20190715_with_blacklist_and_refined"
      
    }
    
    return(l)
    
  })
  
  # Dendrogram tab content ----
  
  # Generate the input dataframe for the bubbleplot
  bubble_input <- reactive({
    
    # Display up to the first 6 genes input
    bubble_prep(gene  = head(input_new()$gene, 6),
                scale = input_new()$scale)
    
  })
  
  # Generate the bubbleplot
  output$bubble <- renderPlot({
    
    req(bubble_input())
    
    bubble_plot(df = bubble_input(),
                max_point_size = input_new()$size)
    
  },
  
  # Choose width to align horizontally with dendrogram image
  width = 1175,
  
  # Customize the height of the bubbleplot to scale with the number of genes which
  # are being displayed, after allocating a baseline height for the x-axis & legend
  height = function() 150 + 30 * length(input_new()$gene))
  
  # Create a tooltip with cluster / expression information that appears when
  # hovering over a bubble
  # 
  # This adapted from this example https://gitlab.com/snippets/16220
  output$bubble_hover_info <- renderUI({
    
    hover <- input$bubble_hover
    
    # Find the nearest data point to the mouse hover position
    point <- nearPoints(bubble_input(),
                        hover,
                        xvar = "Cluster",
                        yvar = "Gene_padded",
                        maxpoints = 1) %>% 
      select(Gene, Sample, Cluster, Cell_type, Cell_class, N_cells,
             Expression, Pct1, Colour)
    
    # Hide the tooltip if mouse is not hovering over a bubble
    if (nrow(point) == 0) return(NULL)
    
    pos <- get_tooltip_pos(hover)
    
    # Create style property fot tooltip
    # background color is set to the cluster colour, with the tooltip a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: ", point$Colour, "cc;",
                    "left:", pos$left_px + 2, "px; top:", pos$top_px + 2, "px; width: 350px;")
    
    # Actual tooltip created as wellPanel, specify info to display
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Gene: </b>",       point$Gene, "<br/>",
                    "<b> Cluster: </b>",    point$Cluster, "<br/>",
                    "<b> Cell type: </b>",  point$Cell_type, "<br/>",
                    "<b> Sample: </b>",     point$Sample, "<br/>",
                    "<b> Expression: </b>", point$Pct1 * point$N_cells, " ",
                    point$Gene, "+ cells out of ", point$N_cells, " cells in cluster <br/>")))
    )
  })
  
  # Download data in bubbleplot tab as TSV
  # output$download_bubble <- downloadHandler(filename = "mean_cluster_expression.tsv",
  #                                           contentType = "text/tsv",
  #                                           content = function(file) {
  #                                             write_tsv(bubble_input() %>% select(-Gene_padded), path = file)
  #                                           })
  
  # Show table with cluster & expression info below bubble plot
  output$cluster_table <- renderDataTable({
    
    req(bubble_input())
    
    bubble_input() %>%
      select(-Pct1, -Gene_padded) %>% 
      mutate(Expression = round(Expression, 2)) %>% 
      spread(Gene, Expression) %>% 
      DT::datatable(options = list(
        columnDefs = list(list(visible = FALSE,
                               # Hide the Colour column
                               targets = c(6))),
        selection = "none")
      ) %>%
      
      # Colour the cluster column based on the palette
      formatStyle("Cluster",
                  backgroundColor = styleEqual(names(joint_mouse_palette), unname(joint_mouse_palette)))
    
  })
  
  # Timecourse tab content ----
  
  # Generate ribbon plot and save the output so that we can later split
  # into the plot itself, and the legend
  ribbon <- reactive({
    
    ribbon_plot(gene   = input_new()$gene[1],
                region = input_new()$region)
    
  })
  
  # Grabbing only the plot part, remove the legend
  output$plotRibbon <- renderPlot({ ribbon() +
      theme(legend.position = "none")
  })
  
  # Extract the ribbon plot legend to plot separately
  output$ribbonLegend <- renderPlot({
    
    leg <- cowplot::get_legend(ribbon())
    plot_grid(leg)
    
  })
  
  # Joint analysis tab content ----
  
  dr_joint_embedding <- reactive({
    
    req(input_new())
    
    region <- input_new()$region
    
    # Load the Cell barcode, 2D coordinates, and selected clustering solution
    df <- feather::read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                                columns = c("Cell",
                                            input_new()$dr,
                                            input_new()$clust))
    
    names(df)[4] <- "Cluster"
    
    return(df)
    
  })
  
  output$dr_joint <- renderPlot({
    
    req(input_new())
    
    dr_plot(dr_joint_embedding(),
            colour_by = "Cluster",
            legend    = FALSE,
            
            # Parameters available to the user
            colours   = input_new()$clust_palette,
            label     = input_new()$label_clusters)
    
  })
  
  output$dr_joint_hover_info <- renderUI({
    
    hover <- input$dr_joint_hover
    
    # Find the nearest data point to the mouse hover position
    point <- nearPoints(dr_joint_embedding(),
                        hover,
                        xvar = input_new()$dr[1],
                        yvar = input_new()$dr[2],
                        maxpoints = 1) %>% 
      select(Cell, Cluster)
    
    # Hide the tooltip if mouse is not hovering over a bubble
    if (nrow(point) == 0) return(NULL)
    
    pos <- get_tooltip_pos(hover)
    
    # Create style property fot tooltip
    # background color is set to the cluster colour, with the tooltip a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", pos$left_px + 2, "px; top:", pos$top_px + 2, "px; width: 350px;")
    
    # Actual tooltip created as wellPanel, specify info to display
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Cell: </b>",    point$Cell, "<br/>",
                    "<b> Cluster: </b>", point$Cluster, "<br/>")))
    )
  })
  
  dr_joint_exp <- reactive({
    
    req(input_new())

    region <- input_new()$region
    
    # Add the gene expression levels to the embedding
    cbind(dr_joint_embedding(),
          feather::read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                                columns = input_new()$gene))

  })
  
  # If more than one gene was provided, compute an aggregate
  dr_joint_agg <- reactive({
    
    req(input_new())
    
    if (length(input_new()$gene) == 1) {
      
      df <- dr_joint_exp()[, 1:5]
      names(df)[5] <- "Expression"
      
    } else if (length(input$gene) > 1) {
      
      # Take the mean of all the gene columns
      meanexp <- dr_joint_exp()[, 5:ncol(dr_joint_exp())] %>% rowMeans
      
      df <- dr_joint_embedding()
      df$Expression <- meanexp
      
    }
    
    return(df)
    
  })
  
  output$feature_joint <- renderPlot({
    
    req(input_new())
    
    feature_plot(dr_joint_agg(),
            label = FALSE,
            
            # Parameters available to the user
            palette = input_new()$ft_palette)
    
  })
  
  output$feature_joint_hover_info <- renderUI({
    
    hover <- input$feature_joint_hover
    
    # Find the nearest data point to the mouse hover position
    point <- nearPoints(dr_joint_agg(),
                        hover,
                        xvar = input_new()$dr[1],
                        yvar = input_new()$dr[2],
                        maxpoints = 1) %>% 
      select(Cell, Cluster, Expression)
    
    # Hide the tooltip if mouse is not hovering over a bubble
    if (nrow(point) == 0) return(NULL)
    
    pos <- get_tooltip_pos(hover)
    
    # Create style property fot tooltip
    # background color is set to the cluster colour, with the tooltip a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", pos$left_px + 2, "px; top:", pos$top_px + 2, "px; width: 350px;")
    
    # Actual tooltip created as wellPanel, specify info to display
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Cell: </b>",    point$Cell, "<br/>",
                    "<b> Cluster: </b>", point$Cluster, "<br/>",
                    "<b> Expression: </b>", round(point$Expression, 2), "<br/>")))
    )
  })
  
  
}

