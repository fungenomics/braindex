
# For the plot cache
# NOTE:
# To get this to work, need to:
# 1) Make the cache folder, i.e. mkdir cache
# 2) On the server, give the shiny user permissions, i.e. chmod -R a=rwx cache
shinyOptions(cache = diskCache("./cache"))

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
      "label_clusters"   = input$label_clusters,
      "ft_palette"       = input$feature_palette,
      "vln_points" = input$vln_points
    )
    
    # Get the columns for the appropriate type of dim red
    if      (input$dr == "tSNE") l$dr <- c("tSNE_1", "tSNE_2")
    else if (input$dr == "UMAP") l$dr <- c("UMAP1", "UMAP2")
    else if (input$dr == "PCA")  l$dr <- c("PC1", "PC2")
    
    # Get the samples depending on region
    if      (input$region == "joint_cortex") l$samples <- c("ct_e12", "ct_e15", "ct_p0", "ct_p3", "ct_p6")
    else if (input$region == "joint_pons")   l$samples <- c("po_e12", "po_e15", "po_p0", "po_p3", "po_p6")
    
    # Get cluster column for each sample depending on region
    if      (input$region == "joint_cortex") l$clust_sample <- "ID_20190730_with_blacklist_and_refined"
    else if (input$region == "joint_pons")   l$clust_sample <- "ID_20190715_with_blacklist_and_refined"
    
    l$clust_palette_sample <- joint_mouse_palette_refined
    
    # Get the clustering to use for the joint analysis tab
    # Option 1: Clustering done at the joint analysis level
    # Option 2: Clustering done per sample (for the forebrain, there was some refinement done so the latest
    # version of labels is different than in the pons)
    if (input$dr_clustering == "joint") {
      
      l$clust <- "ID_20190715_joint_clustering"
      
      if      (input$region == "joint_cortex") l$clust_palette <- cortex_palette_joint
      else if (input$region == "joint_pons")   l$clust_palette <- pons_palette_joint
      
    } else if (input$dr_clustering == "sample") {
      
      l$clust_palette <- joint_mouse_palette_refined
      
      if      (input$region == "joint_cortex") l$clust <- "ID_20190730_with_blacklist_and_refined"
      else if (input$region == "joint_pons")   l$clust <- "ID_20190715_with_blacklist_and_refined"
      
    } else if (input$dr_clustering == "timepoint") {
      
      l$clust <- "orig.ident"
      l$clust_palette <- timepoint_palette
      
    }
    
    return(l)
    
  })
  
  #### ---- Dendrogram tab content ----
  
  # Generate the input dataframe for the bubbleplot 
  bubble_input <- reactive({
    
    # Display up to the first 6 genes input
    # TODO: test 7 - 12 bubble plots. Goal is to display up to 20 together
    bubble_prep(gene  = head(input_new()$gene, 12),
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
  
  # Create a tooltip with cluster / expression information 
  # that appears when hovering over a bubble 
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
    
    # Create style property for tooltip
    # background is set to the cluster colour, with opacity = 100% ("FF" at end of hex)
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: ", point$Colour, "FF;",
                    "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 1, "px; width: 350px;")
    
    # text colour set to white if the background colour is dark, else remains default black
    if (dark(point$Colour)) {
      style <- paste0(style, "color: #FFFFFF")
    }
    
    
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
  
  #### ---- Expression table tab content ----

  # Show table with cluster & expression info 
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
  
  # Download data in bubbleplot tab and expression table as TSV
  output$download_bubble <- downloadHandler(filename = "mean_cluster_expression.tsv",
                                            contentType = "text/tsv",
                                            content = function(file) {
                                              write_tsv(bubble_input() %>% select(-Gene_padded), path = file)
                                            })
  
  #### ---- Timecourse tab content ----
  
  # Generate ribbon plot and save the output so that we can later split
  # into the plot itself, and the legend
  ribbon <- reactive({
    
    ribbon_plot(gene   = input_new()$gene[1],
                region = input_new()$region,
                make_plotly = TRUE)
    
  })
  
  # Add hover functionality to the filled regions of the plot
  output$plotRibbon <- renderPlotly({ 
    add_trace(ribbon(),
              type = "scatter",
              mode = "markers",
              fill = "toself",
              hoveron = "points+fills",
              text = "Points + Fills",
              hoverinfo = "text") 
  })
  
  # Extract the ribbon plot legend to plot separately
  # output$ribbonLegend <- renderPlot({
  #   
  #   leg <- cowplot::get_legend(ribbon())
  #   plot_grid(leg)
  #   ggarrange(leg)
  # })
  
  #### ---- Joint analysis by region tab content ----
  
  dr_joint_embedding <- reactive({
    
    req(input_new())
    
    # Load the Cell barcode, 2D coordinates, and selected clustering solution
    get_embedding(sample  = input_new()$region,
                  dr_cols = input_new()$dr,
                  cluster_column = input_new()$clust)
    
  })
  
  dr_joint <- reactive({
    
    req(input_new())
    
    dr_plot(dr_joint_embedding(),
            colour_by = "Cluster",
            legend    = FALSE,
            
            # Parameters available to the user
            colours   = input_new()$clust_palette,
            label     = input_new()$label_clusters)
    
  })
  
  # output$dr_joint_hover_info <- renderUI({
  #   
  #   hover <- input$dr_joint_hover
  #   
  #   # Find the nearest data point to the mouse hover position
  #   point <- nearPoints(dr_joint_embedding(),
  #                       hover,
  #                       xvar = input_new()$dr[1],
  #                       yvar = input_new()$dr[2],
  #                       maxpoints = 1) %>% 
  #     select(Cell, Cluster)
  #   
  #   # Hide the tooltip if mouse is not hovering over a bubble
  #   if (nrow(point) == 0) return(NULL)
  #   
  #   # Create style property for tooltip
  #   # background color is set to the cluster colour, with the tooltip a bit transparent
  #   # z-index is set so we are sure are tooltip will be on top
  #   style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
  #                   "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 2,  "px; width: 350px;")
  #   
  #   # Actual tooltip created as wellPanel, specify info to display
  #   wellPanel(
  #     style = style,
  #     p(HTML(paste0("<b> Cell: </b>",    point$Cell, "<br/>",
  #                   "<b> ", ifelse(input$dr_clustering == "timepoint", "Sample", "Cluster"),
  #                   "</b>", point$Cluster, "<br/>")))
  #   )
  # })
  
  dr_joint_exp <- reactive({
    
    req(input_new())
    
    get_expression(sample    = input_new()$region,
                   embedding = dr_joint_embedding(),
                   gene      = input_new()$gene,
                   
                   # If more than one gene was provided, compute an aggregate
                   aggregate = TRUE)
    
  })
  
  feature_joint <- reactive({
    
    req(input_new())
    
    feature_plot(dr_joint_exp(),
                 label = FALSE,
                 
                 # Parameters available to the user
                 palette = input_new()$ft_palette)
    
  })
  
  # output$feature_joint_hover_info <- renderUI({
  #   
  #   hover <- input$feature_joint_hover
  #   
  #   # Find the nearest data point to the mouse hover position
  #   point <- nearPoints(dr_joint_agg(),
  #                       hover,
  #                       xvar = input_new()$dr[1],
  #                       yvar = input_new()$dr[2],
  #                       maxpoints = 1) %>% 
  #     select(Cell, Cluster, Expression)
  #   
  #   # Hide the tooltip if mouse is not hovering over a bubble
  #   if (nrow(point) == 0) return(NULL)
  # 
  #   # Create style property fot tooltip
  #   # background color is set to the cluster colour, with the tooltip a bit transparent
  #   # z-index is set so we are sure are tooltip will be on top
  #   style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
  #                   "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 2, "px; width: 350px;")
  #   
  #   # Actual tooltip created as wellPanel, specify info to display
  #   wellPanel(
  #     style = style,
  #     p(HTML(paste0("<b> Cell: </b>",    point$Cell, "<br/>",
  #                   "<b> ", ifelse(input$dr_clustering == "timepoint", "Sample", "Cluster"),
  #                   "</b>", point$Cluster, "<br/>",
  #                   "<b> Expression: </b>", round(point$Expression, 2), "<br/>")))
  #   )
  # })
  
  output$scatter_joint <- renderPlot({
    
    plot_grid(dr_joint(), feature_joint(), rel_widths = c(0.455, 0.545))
    
  })
  
  output$vln_joint <- renderPlot({
    
    vln(dr_joint_exp(),
        palette = input_new()$clust_palette,
        points  = input_new()$vln_points)
    
  })
  
  # ---- Joint analysis by sample tab content ----
  
  dr_sample_embedding <- reactive({
    
    region <- str_split(input_new()$region, "_") %>% sapply(getElement, 2)
    
    df <- cbind(feather::read_feather(glue("data/joint_{region}/{region}.per_sample_embeddings.feather"),
                                      columns = c("Cell", "tSNE_1", "tSNE_2")),
                feather::read_feather(glue("data/joint_{region}/joint_{region}.embedding_and_genes.feather"),
                                      c(input_new()$clust_sample, "orig.ident")))
    
    names(df)[4] <- "Cluster"  
    
    df <- df %>% mutate(Cluster = gsub("_BLACKLISTED", "", Cluster)) %>% 
      tidyr::separate(Cluster, into = c("Prefix", "Cluster"), extra = "drop", sep = "_") %>% 
      select(-Prefix)
    
    dr_sample_list <- split( df, f = df$orig.ident )
    
    return(dr_sample_list)
    
  })
  
  output$dr_sample <- renderPlot({
    
    # Shorten labels for palette
    pal <- input_new()$clust_palette_sample
    names(pal) <- str_split(names(pal), "_") %>% sapply(getElement, 2)
    
    timepoints <- c("E12.5", "E15.5", "P0", "P3", "P6")
    
    map2(dr_sample_embedding(), timepoints,
        ~ dr_plot(.x,
                  colour_by  = "Cluster",
                  legend     = FALSE,
                  hide_ticks = TRUE,
                  title      = .y,
                  label_size = 3,
                  
                  # Parameters available to the user
                  colours   = pal,
                  label     = input_new()$label_clusters,)) %>% 
                  {plot_grid(plotlist = ., ncol = 5)}
    
  })
  
  dr_sample_exp <- reactive({
    
    region <- input_new()$region
    
    df_exp <- cbind(dr_joint_exp() %>% select(Expression),
                    feather::read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                                          c("orig.ident")))
    
    df_sample_exp <- split( df_exp, f = df_exp$orig.ident )
    
    map2(dr_sample_embedding(), df_sample_exp,
         ~ cbind(.x, .y %>% select(- orig.ident)))
    
  })
  
  output$feature_sample <- renderPlot({
    
    map(dr_sample_exp(),
        ~ feature_plot(.x,
                       label = FALSE,
                       hide_ticks = TRUE,
                       
                       # Parameters available to the user
                       palette = input_new()$ft_palette) +
          theme(legend.position = "bottom")) %>%
          {plot_grid(plotlist = ., ncol = 5)}
    
  })
  
  output$vln_sample <- renderPlot({

    # Shorten labels for palette
    pal <- input_new()$clust_palette_sample
    names(pal) <- str_split(names(pal), "_") %>% sapply(getElement, 2)
    
    timepoints <- c("E12.5", "E15.5", "P0", "P3", "P6")
    
    map(dr_sample_exp(),
        ~ vln(.x,
              palette = pal,
              points = input_new()$vln_points) +
          theme(plot.margin = unit(c(0.5, 0, 1, 1.5), "cm"))) %>%
      {plot_grid(plotlist = ., ncol = 1, align = "hv",
                 labels = timepoints, label_size = 18)}
    
  })
  
}

