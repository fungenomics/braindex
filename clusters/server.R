
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
  
  updateSelectizeInput(session, inputId = "gene", choices = genes_anno,
                       server = TRUE)
  
  # Capture all input from this tab as a list in case we want to add
  # more options in the future
  input_new <- eventReactive(input$update, {
    
    g_list <- reactive({
      
      req(input$genelist)

      ext <- tools::file_ext(input$genelist$name)
      switch(ext,
             csv = scan(input$genelist$datapath, 
                        what = "string", sep = ",", 
                        encoding = "UTF-8", fileEncoding = "UTF-8-BOM"),
             tsv = scan(input$genelist$datapath, 
                        what = "string", sep = "\t", 
                        encoding = "UTF-8", fileEncoding = "UTF-8-BOM"),
             txt = scan(input$genelist$datapath, 
                        what = "string", sep = "\t", 
                        encoding = "UTF-8", fileEncoding = "UTF-8-BOM"),
             validate("\n\n\nInvalid file; Please upload a .txt, .csv, or .tsv file")
      )
    })
    
    # Condition which input is used based on the upload toggle
    if (input$upload){
      genes = g_list()
    } else {
      genes = input$gene
    }
    
    # Inputs to use as is
    l <- list(
      "gene"   = genes,
      "scale"  = input$bubble_scale,
      "size"   = input$bubble_size,
      "region" = input$region,
      "label_clusters"   = input$label_clusters,
      "ft_palette"       = input$feature_palette,
      "vln_points" = input$vln_points,
      "plotly_ribbon" = input$plotly_ribbon,
      "mean_exp" = input$mean_exp,
      "heatmap_cells" = input$heatmap_cells,
      "heatmap_anno" = input$heatmap_anno
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
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\n\nPlease enter a gene.")
    )
    
    # Check first 20 gene inputs against the dataset & annotations
    not_data_genes <- check_genes(input_new()$gene, 20, annotation = FALSE)
    not_anno_genes <- check_genes(input_new()$gene, 20, annotation = TRUE)
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }

    validate(
      need(is.null(anno_genes),
           glue("\n\n\n\nThe input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\n\nThe input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Only display mean if more than one gene is given AND the user requested it
    valid_mean <- FALSE
    if (length(input_new()$gene) > 1 && input_new()$mean_exp){
      valid_mean <- TRUE
    } 
    
    # Display the first 20 genes provided as input
    bubble_prep(gene  = head(input_new()$gene, 20),
                scale = input_new()$scale,
                show_mean = valid_mean)
    
  })
  
  # Generate the bubbleplot
  output$bubble <- renderPlot({
    
    req(bubble_input())
    
    bubble_plot(df = bubble_input(),
                max_point_size = input_new()$size)$plot # Get plot part of output
    
  },
  
  # Choose width to align horizontally with dendrogram image
  width = 1103,
  
  # Customize the height of the bubbleplot to scale with the number of genes which
  # are being displayed, after allocating a baseline height for the x-axis & legend
  height = function() 150 + 30 * length(input_new()$gene))
  
  # Create a tooltip with cluster / expression information 
  # that appears when hovering over a bubble 
  # This was adapted from this example: https://gitlab.com/snippets/16220
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
    # background is set to the cluster colour, with opacity = 95% ("F2" at end of hex)
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: ", point$Colour, "F2;",
                    "left: -350px; top: 500px; width: 350px;")
    
    # Set text to white if the background colour is dark, else it's black (default)
    if (dark(point$Colour)) {
      style <- paste0(style, "color: #FFFFFF")
    }
    
    # Specify text content of tooltips - special content for mean expression plot
    if(identical(point$Gene, "MEAN")){
      tooltip_text <- paste0("<b> Mean expression level of plotted genes </b> <br/>",
                             "<b> Cluster: </b>",    point$Cluster, "<br/>",
                             "<b> Cell type: </b>",  point$Cell_type, "<br/>",
                             "<b> Sample: </b>",     point$Sample, "<br/>",
                             "<b> Mean expression: </b>", round(point$Expression, digits=2), "<br/>")
    } else {
      tooltip_text <- paste0("<b> Gene: </b>",       point$Gene, "<br/>",
                             "<b> Cluster: </b>",    point$Cluster, "<br/>",
                             "<b> Cell type: </b>",  point$Cell_type, "<br/>",
                             "<b> Sample: </b>",     point$Sample, "<br/>",
                             "<b> Expression: </b>", point$Pct1 * point$N_cells, " ",
                             point$Gene, "+ cells out of ", point$N_cells, " cells in cluster <br/>")
    }
    
    # Actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(tooltip_text))
    )
  })
  
  # Render the bubble plot gene labels separately with ggdraw
  output$bubble_labels <- renderPlot({
    
    ggdraw(bubble_plot(df = bubble_input(),
                max_point_size = input_new()$size)$labels) # Get labels part of output
    
  },
  
  # Set height of bubble plot gene labels to (hopefully) align with plots
  height = function() 28.5 + 29 * length(input_new()$gene),
  
  # Max length of a gene is 200px
  # NOTE: If altering this, also change the corresponding cellWidth for 
  # splitLayout in ui.R
  width = 200
  
  )
  
  
  #### ---- Expression table tab content ----

  # Show table with cluster & expression info 
  output$cluster_table <- renderReactable({
    req(bubble_input())
    
    # Use the order from bubble_input except reversed 
    gene_table_order <- rev(unique(bubble_input()$Gene))
    
    table <- 
      bubble_input() %>%
      select(-Pct1, -Gene_padded) %>% 
      mutate(Expression = round(Expression, 2)) %>% 
      spread(Gene, Expression) %>% 
      # Select all except Colour column, rename some variables for clarity, and
      # follow bubble_input order for gene columns (saved above)
      select(Cluster,
             Sample,
             "Cell type" = Cell_type,
             "Cell class" = Cell_class,
             "Number of cells" = N_cells,
             all_of(gene_table_order))
    
    # Move mean expression to the rightmost column
    # if ("MEAN" %in% gene_table_order) {
    #   table <- table %>% relocate("MEAN",
    #                              .after = last_col())
    # }
    
    # Create palette for expression level
    orange_pal <- function(x) rgb(colorRamp(c("#ffe4cc", "#ffb54d"))(x), maxColorValue = 255)
    
    # Produce a data table
    reactable(table, 
              rownames = FALSE,
              highlight = TRUE,
              compact = TRUE,
              searchable = TRUE,
              showSortable = TRUE,
              fullWidth = FALSE,
              showPageSizeOptions = TRUE, pageSizeOptions = c(10, 20, 40), defaultPageSize = 10,
              # Formatting for gene columns - color based on expression level
              defaultColDef = colDef(minWidth = 80,
                                     style = function(value) {
                                       color <- orange_pal(value)
                                       list(background = color)
                                     }),
              # Override colDef manually for the first few rows (not genes)
              columns = list(
                Cluster = colDef(minWidth = 110,
                                 style = function(index){
                                   # Colour cluster column background by existing palette
                                   b_color <- toString(unname(joint_mouse_palette)[index])
                                   # Change text colour to white if background is dark
                                   if (dark(b_color)){
                                     f_color = "#FFFFFF"
                                   } else {
                                     f_color = "#000000"
                                   }
                                   list(background = b_color, color = f_color, fontWeight = "bold") 
                                   # # Make the cluster column "sticky" i.e. freeze it in horizontal scroll
                                   #      position = "sticky", left = 0, zIndex = 1)
                                  },
                                 # headerStyle = 
                                 #   list(position = "sticky", left = 0, background = "#fff", zIndex = 1)
                                 ), 
                Sample = colDef(minWidth = 125, style = list(background = "#FFFFFF")),
                "Cell type" = colDef(minWidth = 200, style = list(background = "#FFFFFF")),
                "Cell class" = colDef(minWidth = 150, style = list(background = "#FFFFFF")),
                "Number of cells" = colDef(minWidth = 100, style = list(background = "#FFFFFF"))
                )
    ) 
  })
  
  # output$x4 = renderPrint({
  #   s = input$cluster_table_rows_selected
  #   if (length(s)) {
  #     cat('These clusters were selected:\n\n')
  #     cat(bubble_input()$Cluster[s], sep = ', ')
  #   }
  # })
  
  # Download data in bubbleplot tab and expression table as TSV
  output$download_bubble <- 
    downloadHandler(filename = "mean_cluster_expression.tsv",
                    contentType = "text/tsv",
                    content = function(file) {
                      write_tsv(bubble_input() %>% select(-Gene_padded), path = file)
                    })
  
  #### ---- Timecourse tab content ----
  
  observe({
    x <- input_new()$gene
    
    # Use character(0) to remove all choices when input hasn't been received
    if (is.null(x))
      x <- character(0)
    
    if (length(x) == 1){
      text_select_input <- " input)"
    } else if (length(x) > 1){
      text_select_input <- " inputs)"
    } else {
      text_select_input <- NULL
    }
    
    # Update the label based on length of input
    updateSelectInput(session, "pick_timecourse",
                      label = paste("Select gene to display (from ", length(x), text_select_input),
                      choices = x,
                      selected = head(x, 1) # Select first one by default
    )
  })
  
  # STATIC TIMECOURSE 
  
  # Generate ribbon plot and save the output so that we can allow the
  # user to download it 
  ribbon_static <- reactive({
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\nPlease enter a gene.")
    )
    
    # Check the selected gene against the dataset & annotations
    not_data_genes <- check_genes(input$pick_timecourse, 1, annotation = FALSE)
    not_anno_genes <- check_genes(input$pick_timecourse, 1, annotation = TRUE)
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }
    
    validate(
      need(is.null(anno_genes),
           glue("\n\n\nThe input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\nThe input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Check if expression is all zero in the brain region
    # TODO: fix this to make it work
    all_zero <- ribbon_plot(gene   = input$pick_timecourse,
                      region = input_new()$region)$zero
    
    validate(
      need(all_zero == FALSE, "This gene has no detected expression in the selected brain region.")
    )
    
    p1 <- ribbon_plot(gene   = input$pick_timecourse,
                      region = input_new()$region)$plot
    
    # Get legend using cowplot
    leg <- cowplot::get_legend(p1)
    
    # Grabbing only the plot part, remove the legend
    p1 <- p1 +
    theme(legend.position = "none")
    
    # Combine plot and custom legend into one plot for output
    plot_grid(p1, leg, ncol = 1, rel_heights = c(0.55, 0.45))
    
  })
  
  # Plot leaving some space between x-axis and legend
  output$plotRibbon <- renderPlot({
    
    ribbon_static() +
      theme(plot.margin = unit(c(0, 0, 1, 0), "lines"))
    
  })
  
  # INTERACTIVE TIMECOURSE
  
  # Generate interactive ribbon plot and save the output
  ribbon_plotly <- reactive({

    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\nPlease enter a gene.")
    )
    
    # Check the selected gene against the dataset & annotations
    not_data_genes <- check_genes(input$pick_timecourse, 1, annotation = FALSE)
    not_anno_genes <- check_genes(input$pick_timecourse, 1, annotation = TRUE)
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }
    
    validate(
      need(is.null(anno_genes),
           glue("\n\n\nThe input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\nThe input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Check if expression is all zero in the brain region
    # TODO: fix this to make it work
    all_zero = ribbon_plot(gene   = input$pick_timecourse,
                           region = input_new()$region)$zero
    
    validate(
      need(all_zero == FALSE, "This gene has no detected expression in the selected brain region.")
    )
    
    ribbon_plot(gene   = input$pick_timecourse,
                region = input_new()$region,
                make_plotly = TRUE)$plot
    
  })
  
  output$plotlyRibbon <- renderPlotly({ 

      # Position legend to the right of the plot
      layout(ribbon_plotly(), legend = list(x = 1, y = 0))
  })

  # DOWNLOAD TIMECOURSE (static plot) AS A PDF
  
  output$download_ribbon <- 
    downloadHandler(filename = "timecourse_ribbon.pdf",
                    content = function(file) {
                      ggsave(file, 
                             ribbon_static(),
                             width = 5,
                             height = 5, 
                             units = "in", 
                             scale = 3)
                    },
                    contentType = "application/pdf")
 
  
  #### ---- Region joint analysis tab content ----
  
  dr_joint_embedding <- reactive({
    
    req(input_new())
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\nPlease enter a gene.")
    )
    
    # Check ALL gene inputs against the dataset & annotations
    not_data_genes <- check_genes(input_new()$gene, annotation = FALSE)
    not_anno_genes <- check_genes(input_new()$gene, annotation = TRUE)
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }
    
    validate(
      need(is.null(anno_genes),
           glue("\n\n\nThe input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\nThe input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Load the Cell barcode, 2D coordinates, and selected clustering solution
    get_embedding(sample  = input_new()$region,
                  dr_cols = input_new()$dr,
                  cluster_column = input_new()$clust)
    
  })
  
  # Create a dim. red. plot coloured by cluster
  output$dr_joint <- renderPlot({
    
    req(input_new())
    
    dr_plot(dr_joint_embedding(),
            colour_by = "Cluster",
            legend    = FALSE,
            
            # Parameters available to the user
            colours   = input_new()$clust_palette,
            label     = input_new()$label_clusters) %>% 
      # Get the plot part of list output
      .$plot 
    
  })
  
  clust_centers <- reactive({
    
    req(input_new())
    
    # Get colour associated with each cluster
    palette_df <- tibble::enframe(input_new()$clust_palette, 
                          name = "Cluster", 
                          value = "Colour")
    
    dr_plot(dr_joint_embedding(),
            colour_by = "Cluster",
            legend    = FALSE,
            
            # Parameters available to the user
            colours   = input_new()$clust_palette,
            label     = input_new()$label_clusters) %>% 
      # Get the centers part of list output
      .$centers %>% 
      # Add cluster colours into the same dataframe
      left_join(palette_df, by = c("Cluster" = "Cluster")) 
      
    
  })
  
  output$dr_joint_hover_info <- renderUI({

    hover <- input$dr_joint_hover

    # Find the nearest data point to the mouse hover position
    point <- nearPoints(clust_centers(),
                        hover,
                        xvar = "center_x",
                        yvar = "center_y",
                        threshold = 25,
                        maxpoints = 1) %>%
      select(Cluster, Colour)

    # Hide the tooltip if the cluster information is blank
    if (nrow(point) == 0) return(NULL)

    # Create style property for tooltip
    # background is set to the cluster colour, with opacity = 95% ("F2" at end of hex)
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color:", point$Colour, "F2;",
                    "left:", hover$coords_css$x, "px; top:", hover$coords_css$y + 175,  "px; 
                    width: auto; height: 60px;")
    
    # Set text to white if the background colour is dark, else it's black (default)
    if (dark(point$Colour)) {
      style <- paste0(style, "color: #FFFFFF")
    }
    
    # Actual tooltip created as wellPanel, specify info to display
    wellPanel(
      style = style,
      p(HTML(paste0(
        #"<b> Cell: </b>",    point$Cell, "<br/>",
                    #"<b> ", ifelse(input$dr_clustering == "timepoint", "Sample", "Cluster"),
                    "<b>Cluster: </b>", point$Cluster)))
    )
  })
  
  # Get the gene expression for the gene(s) from user input
  dr_joint_exp <- reactive({
    
    req(input_new())
    
    express <- get_expression(sample    = input_new()$region,
                              embedding = dr_joint_embedding(),
                              gene      = input_new()$gene,
                              # If more than one gene was provided, compute an aggregate
                              aggregate = TRUE)
    
    # Display message to the user if there is 0 expression throughout region
    validate(
      need(express$zero == FALSE, "This gene has no detected expression in the selected brain region.")
    )
    
    express$data
    
  })
  
  # Create a dim. red. plot coloured by expression
  output$feature_joint <- renderPlot({
    
    req(input_new())
    
    feature_plot(dr_joint_exp(),
                 label = FALSE,
                 
                 # Parameter available to the user
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
  
  # # Plot the dim. red. plots together
  # output$scatter_joint <- renderPlot({
  #   
  #   plot_grid(dr_joint(), feature_joint(), rel_widths = c(0.455, 0.545))
  #   
  # })
  
  # Create and plot a violin plot coloured by cluster
  output$vln_joint <- renderPlot({
    
    vln(dr_joint_exp(),
        palette = input_new()$clust_palette,
        points  = input_new()$vln_points)
    
  })
  
  #### ---- Sample joint analysis tab content ----
  
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
  
  # Create and plot dim. red. plots for each timepoint, coloured by cluster
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
                  label     = input_new()$label_clusters,
                  
                  # Get the plot part of list output
                  )$plot) %>% 
            {plot_grid(plotlist = ., ncol = 5)}
    
  })
  
  # Get the gene expression for the gene(s) from user input
  dr_sample_exp <- reactive({
    
    region <- input_new()$region
    
    df_exp <- cbind(dr_joint_exp() %>% select(Expression),
                    feather::read_feather(glue("data/{region}/{region}.embedding_and_genes.feather"),
                                          c("orig.ident")))
    
    df_sample_exp <- split( df_exp, f = df_exp$orig.ident )
    
    map2(dr_sample_embedding(), df_sample_exp,
         ~ cbind(.x, .y %>% select(- orig.ident)))
    
  })
  
  # Create and plot dim. red. plots for each timepoint, coloured by expression
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
  
  # Create and plot violin plots for each timepoint, coloured by cluster
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
                 labels = timepoints, label_size = 15)}
    
  })
  
  #### ---- Clusters ranked by expression tab content ----
  
  observe({
    x <- input_new()$gene
    
    # Use character(0) to remove all choices when input hasn't been received
    if (is.null(x))
      x <- character(0)
    
    if (length(x) == 1){
      text_select_input <- " input)"
    } else if (length(x) > 1){
      text_select_input <- " inputs)"
    } else {
      text_select_input <- NULL
    }
    
    # Update the label based on length of input
    updateSelectInput(session, "pick_ranked",
                      label = paste("Select gene to display (from ", length(x), text_select_input),
                      choices = x,
                      selected = head(x, 1) # Select first one by default
    )
  })
  
  ranked_plot <- reactive({
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\n      Please enter a gene.")
    )
    
    if(input_new()$mean_exp){
      # Check ALL inputs against the dataset genes
      num_genes <- NULL
      # Check gene inputs against the dataset & annotations
      not_data_genes <- check_genes(input_new()$gene, num_genes, annotation = FALSE)
      not_anno_genes <- check_genes(input_new()$gene, num_genes, annotation = TRUE)
    } else{
      # Check selected gene input against the dataset & annotations
      not_data_genes <- check_genes(input$pick_ranked, annotation = FALSE)
      not_anno_genes <- check_genes(input$pick_ranked, annotation = TRUE)
    }
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }
    
    validate(
      need(is.null(anno_genes),
           glue("\n\n\n     The input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\n     The input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    if (input_new()$mean_exp){
      df <- bubble_prep(gene = input_new()$gene,
                        show_mean = TRUE) %>% 
        filter(Gene == "MEAN")
      y_axis_text <- "Mean gene expression"
      title_text <- "Mean expression over all selected genes"
    } else {
      df <- bubble_prep(gene = input$pick_ranked)
      y_axis_text <- glue("Mean {input$pick_ranked} expression")
      title_text <- input$pick_ranked
    }
    
    df <- df %>% 
      # Order from highest to lowest by expression (ranked)
      arrange(desc(Expression)) %>% 
      mutate(Cluster = factor(Cluster, levels = .$Cluster)) %>% 
      # Rename cell classes to more general names
      mutate(Cell_class = case_when(
        grepl("RGC", Cell_class) | grepl("-P$", Cluster) ~ "Progenitors/cyc.",
        grepl("Olig", Cell_class) ~ "Oligodendrocytes",
        grepl("Epen", Cell_class) ~ "Ependymal",
        grepl("Astr", Cell_class) ~ "Astrocytes",
        grepl("[Nn]euron", Cell_class) ~ "Neurons",
        grepl("Non-neuro|Immune", Cell_class) ~ "Non-neuroect.",
        TRUE ~ "Other"
      ))
    
    p1 <- df %>% ggplot(aes(x = Cluster, y = Expression)) +
      geom_bar(aes(fill = Cluster), stat = "identity") +
      scale_fill_manual(values = df$Colour) +
      theme_min(border_colour = "gray90") +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            # Remove white space at the bottom of plot
            plot.margin = margin(b=0, unit="cm")) +
      expand_limits(x = -18) +
      labs(title = title_text) +
      ylab(y_axis_text)

    ticks <- ggplot() + add_class_ticks(df, unique(df$Cell_class),
                             palette = general_palette,
                             start = -5, sep = 5, height = 30, label_x_pos = -16, fontsize = 3) +
      # Make sure to expand to the same value that's in p1
      expand_limits(x = -18) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # Remove plot border
            panel.border = element_blank(),
            # Remove white space at the top of the plot
            plot.margin = margin(t=0, unit="cm")) 
    
    plot_grid(p1, ticks, ncol = 1, align = "v")
  })

  output$rank_tick_plot <- renderPlot({
    ranked_plot()
  })
  
  output$download_ranked_plot <- 
    downloadHandler(filename = "ranked_plot.pdf",
                    content = function(file) {
                      ggsave(file, 
                             ranked_plot(),
                             width = 18,
                             height = 6,
                             units = "in", 
                             scale = 1)
                    },
                    contentType = "application/pdf")
  
  #### ---- Cell types clustered by expression tab content ----
  
  # Initial values, will show up briefly while plot is loading, so set 
  # large height value so the plot that shows briefly is out of view
  heatmap_dims <- reactiveValues(height = "100in", width = "12in")
  
  heatmap_plot <- reactive({
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\n    Please enter a gene.")
    )
    
    # Check if number of genes > 1 - need at least 2 genes to cluster
    validate(
      need(length(input_new()$gene) > 1, "\n\n\n    Please enter more than one gene. Heatmap clustering requires at least two genes.")
    )
    
    # Check ALL gene inputs against the dataset & annotations
    not_data_genes <- check_genes(input_new()$gene, annotation = FALSE)
    not_anno_genes <- check_genes(input_new()$gene, annotation = TRUE)
    
    # Store genes that are within the annotation but not in dataset 
    if (setequal(not_data_genes, not_anno_genes)){
      anno_genes <- NULL
    } else {
      anno_genes <- setdiff(not_data_genes, not_anno_genes)
    }
    
    validate(
      need(is.null(anno_genes),
           glue("\n\n\n    The input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\n    The input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Check that at least one cell type has been chosen by the user
    validate(
      need(length(input_new()$heatmap_cells) > 0,
           glue("\n\n\n    Please select at least one cell type."))
    )
    
    # Check that at least one column annotation has been chosen by the user
    validate(
      need(length(input_new()$heatmap_anno) > 0,
           glue("\n\n\n    Please select at least one column annotation."))
    )
    
    # Get gene expression values for input genes
    df <- bubble_prep(gene = input_new()$gene) 
    
    df <- df %>% 
      # Rename cell classes to more general names
      mutate(Cell_class = case_when(
        grepl("RGC", Cell_class) | grepl("-P$", Cluster) ~ "Progenitors/cyc.",
        grepl("Olig", Cell_class) ~ "Oligodendrocytes",
        grepl("Epen", Cell_class) ~ "Ependymal",
        grepl("Astr", Cell_class) ~ "Astrocytes",
        grepl("[Nn]euron", Cell_class) ~ "Neurons",
        grepl("Non-neuro|Immune", Cell_class) ~ "Non-neuroect.",
        TRUE ~ "Other"
      ))
    
    # Add new columns to annotation df for brain region and time points information
    df <- df %>% 
      mutate(Region = case_when(
        grepl("Forebrain", Sample) ~ "Forebrain",
        grepl("Pons", Sample) ~ "Pons",
        TRUE ~ "Other"
      ))
    df <- df %>% 
      mutate(Timepoint = case_when(
        grepl("E12.5", Sample) ~ "E12.5",
        grepl("E15.5", Sample) ~ "E15.5",
        grepl("P0", Sample) ~ "P0",
        grepl("P3", Sample) ~ "P3",
        grepl("P6", Sample) ~ "P6",
        TRUE ~ "Other"
      ))
    
    # Store the df for later use in annotations
    df_for_anno <- df 
    
    # Store mean expression for each cluster 
    df <- df %>% 
      group_by(Gene, Cluster, Cell_class) %>% 
      summarize(Expression = mean(Expression))
    
    df <- df %>% filter(Cell_class %in% input_new()$heatmap_cells)
    
    # Select only certain clusters that can be categorized 
    # (from Selin's code, to be confirmed)
    # df <- df %>% filter(grepl("ASTR|EPEN|OL|OPC|EXN", Cluster) & !grepl("^B", Cluster)) 
    # df <- df %>% filter(!grepl("^B", Cluster)) 
    
    # Pivot the cluster rows to columns
    df <- df %>% 
      select(Gene, Cluster, Expression) %>% 
      mutate(Expression = as.numeric(Expression)) %>% 
      tidyr::pivot_wider(names_from = "Cluster", values_from = "Expression")  
    
    # Set NA values in df to 0 (from Selin's code)
    df[is.na(df)] <- 0 
    
    # Flip dataframe over to match the order of user input (bubble_prep outputs reverse)
    df <- df[nrow(df):1,]
    
    # Convert dataframe values to matrix, set rownames to genes for labeling purposes
    mat <- df[,-1] %>%
      data.matrix()  
    rownames(mat) <- df$Gene
    
    # Create values for heatmap annotations (coloured bars at top)
    
      # CELL CLASS ANNOTATIONS
    hm_anno_class <- makePheatmapAnno(general_palette, "Cell_class")
    hm_anno_class$anno_row <- left_join(hm_anno_class$anno_row, 
                                  unique(select(df_for_anno, Cluster, Cell_class)), by = "Cell_class") 
    rownames(hm_anno_class$anno_row) <- hm_anno_class$anno_row$Cluster
    hm_anno_class$anno_row$Cluster <- NULL # Prevent individual clusters from showing in plot
      
      # REGION ANNOTATIONS
    hm_anno_region <- makePheatmapAnno(region_palette, "Region")
    hm_anno_region$anno_row <- left_join(hm_anno_region$anno_row,
                                         unique(select(df_for_anno, Cluster, Region), by = "Region"))
    rownames(hm_anno_region$anno_row) <- hm_anno_region$anno_row$Cluster
    hm_anno_region$anno_row$Cluster <- NULL # Prevent individual clusters from showing in plot
      
      # TIMEPOINTS ANNOTATIONS
    hm_anno_time <- makePheatmapAnno(timepoint_palette, "Timepoint")
    hm_anno_time$anno_row <- left_join(hm_anno_time$anno_row,
                                         unique(select(df_for_anno, Cluster, Timepoint), by = "Timepoint"))
    rownames(hm_anno_time$anno_row) <- hm_anno_time$anno_row$Cluster
    hm_anno_time$anno_row$Cluster <- NULL # Prevent individual clusters from showing in plo
    
    anno = data.frame(Cell_class = hm_anno_class$anno_row,
                       Region = hm_anno_region$anno_row,
                       Timepoint = hm_anno_time$anno_row)
    # anno_colors = list(Cell_class = hm_anno_class$side_colors,
    #                    Region = hm_anno_region$side_colors,
    #                    Timepoint = hm_anno_time$side_colors)
    
    # Display only the annotations selected by the user
    anno <- anno %>% select(input_new()$heatmap_anno)
    # anno_colors <- anno_colors %>% select(input_new()$heatmap_anno)
    
    # Plot heatmap, store dimensions for dynamically plotting full width
    hm <- mat %>% 
          apply(1, scales::rescale) %>% 
          t() %>% 
          pheatmap::pheatmap(border_color = NA,
                             color = colorRampPalette(c("blue", "white", "red"))(100),
                             scale = "none",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             cellwidth = 17,
                             cellheight = 17,
                             fontsize = 13,
                             annotation_col = anno, 
                             #annotation_colors = anno_colors, 
                             na_col = "#e5e5e5")

    heatmap_dims$width <- glue("{get_plot_dims(hm)$width}in")
    heatmap_dims$height <- glue("{get_plot_dims(hm)$height}in")
    hm
  })
  
  output$heatmap <- renderPlot(heatmap_plot())
  
  # Plot the heatmap using dynamic width and height values stored above
  output$heatmapUI <- renderUI({
    plotOutput("heatmap", 
               width = heatmap_dims$width, 
               height = heatmap_dims$height)
  })
  
  output$download_heatmap <- 
    downloadHandler(filename = "heatmap.pdf",
                    content = function(file) {
                      ggsave(file, 
                             heatmap_plot(),
                             width = (as.numeric(substr(heatmap_dims$width,
                                                       1, nchar(heatmap_dims$width)-2)) + 1.5),
                             height = (as.numeric(substr(heatmap_dims$height,
                                                        1, nchar(heatmap_dims$height))) + 2.5),
                             units = "in", 
                             scale = 1)
                    },
                    contentType = "application/pdf")
  
}

