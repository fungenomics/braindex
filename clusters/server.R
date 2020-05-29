
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
  
  # Dendrogram tab content ----
  
  # Capture all input from this tab as a list in case we want to add
  # more options in the future
  input_new <- eventReactive(input$update, {
    
    # Inputs to use as is
    l <- list(
      "gene"   = input$gene,
      "scale"  = input$bubble_scale,
      "size"   = input$bubble_size,
      "region" = input$region,
      "label_clusters" = input$label_clusters
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
             Expression, Pct1)
    
    # Hide the tooltip if mouse is not hovering over a bubble
    if (nrow(point) == 0) return(NULL)
    
    # Create tooltip for mouseover
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # Create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px; width: 350px;")
    
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
    
    region <- input_new()$region
    
    df <- feather::read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                                columns = c("Cell",
                                            input_new()$dr,
                                            input_new()$clust))
    
    names(df)[4] <- "Cluster"
    
    return(df)
    
  })
  
  output$dr_joint <- renderPlot({
    
    dr_plot(dr_joint_embedding(),
            colour_by = "Cluster",
            colours   = input_new()$clust_palette,
            legend    = FALSE,
            label     = input_new()$label_clusters)
    
  })
  
  # dr_joint_exp <- reactive({
  #   
  #   req(dr_k27m())
  #   cbind(dr_k27m(),
  #         feather::read_feather("data/k27m.embedding_and_genes.feather", columns = input$gene))
  #   
  # })
  
  
}

