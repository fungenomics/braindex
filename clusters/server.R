
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
  
  #### ---- Sidebar input ----
  
  # Implement selectize (i.e. server-side selection) to speed up loading time
  updateSelectizeInput(session, inputId = "gene", choices = genes_anno,
                       server = TRUE)
  
  # Capture all input from this tab as a list in case we want to add
  # more options in the future
  input_new <- eventReactive(input$update, {
    
    # Allow user to upload gene lists as .csv, .tsv, .txt files
    #   encoded in UTF-8
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
             # Display message if the file is not in one of these formats
             validate("\n\n\nInvalid file; Please upload a .txt, .csv, or .tsv file")
      )
    })
    
    # Switch between gene input by file or typing based on toggle button
    if (input$upload)      genes = g_list()
    else                   genes = input$gene
    
    
    # Use the following inputs without modification:
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
    
    # Modify other inputs as required:
    
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
    # Option 2: Clustering done per sample (for the forebrain, there was some 
    #   refinement done, so the latest version of labels is different than in pons)
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
    
    # Perform input validation
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
    
    # Custom plotter (functions.R)
    bubble_plot(df = bubble_input(),
                max_point_size = input_new()$size)$plot # Get plot part of output
    
  },
  
  # Choose specific width to align horizontally with dendrogram image
  width = 1103,
  
  # Customize height of the bubbleplot to scale with the number of displayed 
  # genes, after allocating a baseline height for the x-axis & legend
  height = function() 150 + 30 * length(input_new()$gene))
  
  # Create a tooltip with cluster / expression information,
  #     which appears when hovering over a bubble 
  # Adapted from this example: https://gitlab.com/snippets/16220 
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
    # background is set to the cluster colour, opacity = 95% ("F2" at end of hex)
    # z-index is set so we are sure the tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: ", point$Colour, "F2;",
                    "left: -350px; top: 500px; width: 350px;")
    
    # Set text to white (override default black) if background colour is dark
    # See dark() in functions.R
    if (dark(point$Colour)) {
      style <- paste0(style, "color: #FFFFFF")
    }
    
    # Specify text content of tooltips, with special content for mean expression plot
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
    
    # Created tooltip as wellPanel
    wellPanel(
      style = style,
      p(HTML(tooltip_text))
    )
  })
  
  # Render bubble plot gene labels separately with ggdraw
  output$bubble_labels <- renderPlot({
    
    ggdraw(bubble_plot(df = bubble_input(),
                max_point_size = input_new()$size)$labels) # Get labels part of output
    
  },
  
  # Set height of bubble plot gene labels to align with plots
  height = function() 28.5 + 29 * length(input_new()$gene),
  
  # Max length of a gene is 200px
  # NOTE: If altering this, also change the corresponding cellWidth for 
  # splitLayout in ui.R
  width = 200
  
  )
  
  
  #### ---- Cluster info & markers table tab content ----

  # --- EXPRESSION TABLE ---
  
  # Display table before update button has been clicked 
  output$cluster_table_no_update <- renderReactable({
      
    # Modify metadata for plotting
      table <- 
        metadata %>%
        as.data.frame() %>%
        select(Cluster,
               Sample,
               "Cell type" = Cell_type,
               "Cell class" = Cell_class,
               "Number of cells" = N_cells)
      
      # Display the dataframe as a reacTable with custom partial function (functions.R)
      reactable_table(table, 
                      # Column formatting
                      columns = 
                        list(Cluster = colDef(minWidth = 110,
                                              style = function(value){
                                                # Colour cluster column background by existing palette
                                                b_color <- toString(filter(metadata, Cluster == value)$Colour)
                                                # Change text colour to white if background is dark
                                                if (dark(b_color)){
                                                  f_color = "#FFFFFF"
                                                } else {
                                                  f_color = "#000000"
                                                }
                                                list(background = b_color, color = f_color, fontWeight = "bold") 
                                                ## Make the cluster column "sticky" i.e. freeze it in horizontal scroll
                                                #position = "sticky", left = 0, zIndex = 1)
                                              },
                                              # headerStyle = 
                                              #   list(position = "sticky", left = 0, background = "#fff", zIndex = 1)
                        ), 
                        
                        Sample = colDef(minWidth = 125),
                        "Cell type" = colDef(minWidth = 200),
                        "Cell class" = colDef(minWidth = 150),
                        "Number of cells" = colDef(minWidth = 100)
                        ),
                      # Implement row selection and formatting
                      selection = "single", onClick = "select", 
                      theme = reactableTheme(
                        rowSelectedStyle = 
                          list(backgroundColor = "#ccc", 
                               boxShadow = "inset 5px 0 0 0 #ffa62d")
                      )
      ) 
  })
  
  # Show table with cluster & expression info 
  output$cluster_table <- renderReactable({
    
    if (length(input_new()$gene) > 0){ 
      
      req(bubble_input()) # Use the reactive value stored when displaying the dendrogram tab
      
      # Store reversed bubble_input gene order (to match user input left to right)
      gene_table_order <- rev(unique(bubble_input()$Gene))
      
      # Modify table for plotting
      table <- 
        bubble_input() %>%
        select(-Pct1, -Gene_padded) %>% 
        mutate(Expression = round(Expression, 2)) %>% 
        spread(Gene, Expression) %>% 
        # Select all except Colour column, rename columns for human readability,
        # follow bubble_input order for gene columns (saved above)
        select(Cluster,
               Sample,
               "Cell type" = Cell_type,
               "Cell class" = Cell_class,
               "Number of cells" = N_cells,
               all_of(gene_table_order))
      
      ## Move mean expression to the rightmost column
      # if ("MEAN" %in% gene_table_order) {
      #   table <- table %>% relocate("MEAN",
      #                              .after = last_col())
      # }
      
      # Store a dynamic list of gene columns (based on input) with assigned formatting
      # Concept from: https://github.com/glin/reactable/issues/65#issuecomment-667577253 
      gene_columns <- 
        sapply(gene_table_order, # Vector of input genes
               function(this_gene) {
                 return(list(this_gene = 
                               colDef(minWidth = 80,
                                      # Color gene column values based on expression level
                                      style = function(value) {
                                        color <- orange_pal(value)
                                        list(background = color)
                                      })))
               })
      
      # Rename list names to gene names (they were changed during sapply)
      names(gene_columns) <- gene_table_order 
      
      # Display the dataframe as a reacTable with custom partial function (functions.R)
      reactable_table(table, 
                      # Column formatting
                      columns = c(
                        list(Cluster = colDef(minWidth = 110,
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
                                           ## Make the cluster column "sticky" i.e. freeze it in horizontal scroll
                                           #position = "sticky", left = 0, zIndex = 1)
                                          },
                                         # headerStyle = 
                                         #   list(position = "sticky", left = 0, background = "#fff", zIndex = 1)
                                         )), 
                        
                        list(Sample = colDef(minWidth = 125)),
                        list("Cell type" = colDef(minWidth = 200)),
                        list("Cell class" = colDef(minWidth = 150)),
                        list("Number of cells" = colDef(minWidth = 100)),
                        
                        # Include dynamic list of columns stored above
                        gene_columns
                      ),
                      # Implement row selection and formatting
                      selection = "single", onClick = "select", 
                      theme = reactableTheme(
                        rowSelectedStyle = 
                          list(backgroundColor = "#ccc", 
                               boxShadow = "inset 5px 0 0 0 #ffa62d")
                      )
      ) 
    
    } else { # DISPLAY CLUSTER TABLE WITHOUT GENES HAVING BEEN ENTERED, AFTER UPDATE BUTTON CLICKED
      
      # Modify metadata for plotting
      table <- 
        metadata %>%
        as.data.frame() %>%
        select(Cluster,
               Sample,
               "Cell type" = Cell_type,
               "Cell class" = Cell_class,
               "Number of cells" = N_cells)
      
      # Display the dataframe as a reacTable with custom partial function (functions.R)
      reactable_table(table, 
                      # Column formatting
                      columns = 
                        list(Cluster = colDef(minWidth = 110,
                                              style = function(value){
                                                # Colour cluster column background by existing palette
                                                b_color <- toString(filter(metadata, Cluster == value)$Colour)
                                                # Change text colour to white if background is dark
                                                if (dark(b_color)){
                                                  f_color = "#FFFFFF"
                                                } else {
                                                  f_color = "#000000"
                                                }
                                                list(background = b_color, color = f_color, fontWeight = "bold") 
                                                ## Make the cluster column "sticky" i.e. freeze it in horizontal scroll
                                                #position = "sticky", left = 0, zIndex = 1)
                                              },
                                              # headerStyle = 
                                              #   list(position = "sticky", left = 0, background = "#fff", zIndex = 1)
                        ), 
                        
                        Sample = colDef(minWidth = 125),
                        "Cell type" = colDef(minWidth = 200),
                        "Cell class" = colDef(minWidth = 150),
                        "Number of cells" = colDef(minWidth = 100)
                      ),
                      # Implement row selection and formatting
                      selection = "single", onClick = "select", 
                      theme = reactableTheme(
                        rowSelectedStyle = 
                          list(backgroundColor = "#ccc", 
                               boxShadow = "inset 5px 0 0 0 #ffa62d")
                      )
      ) 
      
    }
  })
  
  # Download data in bubbleplot tab and expression table as TSV
   
    #reactive({
    # Download expression table if gene(s) have been entered
    #
      output$download_bubble <- 
        # reactive({ 
          # if (length(input_new()$gene) > 0) { 
          downloadHandler(filename = "mean_cluster_expression.tsv",
                          contentType = "text/tsv",
                          content = function(file) {
                            write_tsv(bubble_input() %>% select(-Gene_padded), path = file)
                          })
      #   } else {
      #     downloadHandler(filename = "cluster_information.tsv",
      #                     contentType = "text/tsv",
      #                     content = function(file) {
      #                       write_tsv(as.data.frame(metadata), path = file)
      #                     })
      #   }
      # })
    # } else { # Else just download the metadata (cluster info without genes)
    #   output$download_bubble <- downloadHandler(filename = "cluster_information.tsv",
    #                   contentType = "text/tsv",
    #                   content = function(file) {
    #                     write_tsv(metadata %>% select(-Gene_padded), path = file)
    #                   })
    # }
    # })

  
  # --- MARKER TABLE ---
  
  # Store selected row (cluster) as table index, extract cluster name  
  selected_index <- reactive(getReactableState("cluster_table", "selected"))
  selected_index_no_update <- reactive(getReactableState("cluster_table_no_update", "selected"))
  cluster_order <- read_feather("data/joint_mouse/mean_expression_per_ID_20190715_cluster.feather",
                                columns = c("Cluster")) 
  selected_cluster <- reactive(cluster_order$Cluster[selected_index()])
  selected_cluster_no_update <- reactive(cluster_order$Cluster[selected_index_no_update()])
  
  # Display the selected cluster's name to the user above marker table
  output$selected_clust <- renderUI({
    HTML(glue("<h4>Selected cluster: {selected_cluster()}</h4>"))
  })
  output$selected_clust_no_update <- renderUI({
    HTML(glue("<h4>Selected cluster: {selected_cluster_no_update()}</h4>"))
  })
  
  # Extract selected cluster's info from metadata
  # Beware: metadata also includes some human samples
  cluster_info <- reactive(
    metadata %>% 
      filter(Cluster_nounderscore == as.character(selected_cluster())) %>% 
      select(Alias, Structure, Age, Cluster_number)
  )
  
  cluster_info_no_update <- reactive(
    metadata %>% 
      filter(Cluster_nounderscore == as.character(selected_cluster_no_update())) %>% 
      select(Alias, Structure, Age, Cluster_number)
  )
  
  output$marker_table_no_update <- renderReactable({
    
    # Validate whether a cluster is selected or not, display message to user
    validate(
      need(!is.null(getReactableState("cluster_table_no_update", "selected")), 
             "Please select a cluster for which to display markers in the expression table above.")
    )
    
    # Load marker and signature files for the selected cluster's region & timepoint
    if (as.character(cluster_info_no_update()$Structure) == "Forebrain") region <- "cortex"
    else region <- "pons"  # else includes samples labeled "Pons" as well as "Hindbrain"
    
    markers <- data.table::fread(glue("data/markers_{region}/{cluster_info_no_update()$Alias}.markers.tsv.gz"),
                                 data.table = FALSE)
    signature <- loadRData(glue("data/signatures_{region}/{cluster_info_no_update()$Alias}.signatures_no_mito.Rda"))
    
    # Store selected cluster's number as a string for filtering
    cluster_string <- as.character(cluster_info_no_update()$Cluster_num)
    
    markers_filter <- markers %>% 
      # Filter markers to contain only genes in the specified cluster's signature
      filter(cluster == as.integer(cluster_info_no_update()$Cluster_number)) %>% 
      filter(external_gene_name %in% signature[[2]][[cluster_string]]) %>% 
      
      # Calculate new columns and round decimal values
      mutate(Specificity = pct.1 - pct.2) %>% 
      mutate(avg_logFC = round(avg_logFC, 4)) %>%
      mutate(p_val_adj = signif(p_val_adj, 4)) %>%
      mutate(pct.1 = round(pct.1, 4)) %>%
      mutate(pct.2 = round(pct.2, 4)) %>%
      mutate(Specificity = round(Specificity, 4)) %>%
      
      # Rename columns for human readability
      select(Gene = external_gene_name,
             "Gene type" = gene_biotype,
             Description = description,
             "Average log fold change" = avg_logFC,
             "P-value" = p_val_adj,
             "Detection rate in cluster" = pct.1,
             "Detection rate outside cluster" = pct.2,
             Specificity)
    
    # Output the joined dataframe as a reacTable
    reactable_table(markers_filter,
                    # Column formatting
                    columns = c(
                      list(Gene = colDef(minWidth = 90, style = list(fontWeight = "bold")),
                           "Gene type" = colDef(minWidth = 120
                                                # ,
                                                ## Conditional colour formatting
                                                # style = function(value) {
                                                #   if (is.na(value)) color <- "e5e5e5"
                                                #   else if (value == "protein_coding") color <- "#ebbbab"
                                                #   else color <- "#abd3eb"
                                                #   list(background = color)
                                                # }
                           ),
                           Description = colDef(minWidth = 300),
                           "Average log fold change" = colDef(minWidth = 120),
                           "P-value" = colDef(minWidth = 110),
                           "Detection rate in cluster" = colDef(minWidth = 120),
                           "Detection rate outside cluster" = colDef(minWidth = 130),
                           Specificity = colDef(minWidth = 110,
                                                # Conditional colour formatting
                                                style = function(value) {
                                                  color <- orange_pal(value)
                                                  list(background = color)
                                                })
                      )
  
    ))
  })
  
  output$marker_table <- renderReactable({
    
    # Validate whether a cluster is selected or not, display message to user
    validate(
      need(!is.null(getReactableState("cluster_table", "selected")), "
           Please select a cluster for which to display markers in the expression table above.")
    )
    
    # Load marker and signature files for the selected cluster's region & timepoint
    if (as.character(cluster_info()$Structure) == "Forebrain") region <- "cortex"
    else region <- "pons"  # else includes samples labeled "Pons" as well as "Hindbrain"
    
    markers <- data.table::fread(glue("data/markers_{region}/{cluster_info()$Alias}.markers.tsv.gz"),
                                 data.table = FALSE)
    signature <- loadRData(glue("data/signatures_{region}/{cluster_info()$Alias}.signatures_no_mito.Rda"))
    
    # Store selected cluster's number as a string for filtering
    cluster_string <- as.character(cluster_info()$Cluster_num)
    
    markers_filter <- markers %>% 
      # Filter markers to contain only genes in the specified cluster's signature
      filter(cluster == as.integer(cluster_info()$Cluster_number)) %>% 
      filter(external_gene_name %in% signature[[2]][[cluster_string]]) %>% 
      
      # Calculate new columns and round decimal values
      mutate(Specificity = pct.1 - pct.2) %>% 
      mutate(avg_logFC = round(avg_logFC, 4)) %>%
      mutate(p_val_adj = signif(p_val_adj, 4)) %>%
      mutate(pct.1 = round(pct.1, 4)) %>%
      mutate(pct.2 = round(pct.2, 4)) %>%
      mutate(Specificity = round(Specificity, 4)) %>%
      
      # Rename columns for human readability
      select(Gene = external_gene_name,
             "Gene type" = gene_biotype,
             Description = description,
             "Average log fold change" = avg_logFC,
             "P-value" = p_val_adj,
             "Detection rate in cluster" = pct.1,
             "Detection rate outside cluster" = pct.2,
             Specificity)
    
    # Output the joined dataframe as a reacTable
    reactable_table(markers_filter,
                    # Column formatting
                    columns = c(
                      list(Gene = colDef(minWidth = 90, style = list(fontWeight = "bold")),
                           "Gene type" = colDef(minWidth = 120
                                                # ,
                                                ## Conditional colour formatting
                                                # style = function(value) {
                                                #   if (is.na(value)) color <- "e5e5e5"
                                                #   else if (value == "protein_coding") color <- "#ebbbab"
                                                #   else color <- "#abd3eb"
                                                #   list(background = color)
                                                # }
                                                ),
                           Description = colDef(minWidth = 300),
                           "Average log fold change" = colDef(minWidth = 120),
                           "P-value" = colDef(minWidth = 110),
                           "Detection rate in cluster" = colDef(minWidth = 120),
                           "Detection rate outside cluster" = colDef(minWidth = 130),
                           Specificity = colDef(minWidth = 110,
                                                # Conditional colour formatting
                                                style = function(value) {
                                                  color <- orange_pal(value)
                                                  list(background = color)
                                                })
                    )
    ))
  })
  
  #### ---- Timecourse tab content ----
  
  # Dynamic dropdown menu to choose which input gene to plot
  observe({
    x <- input_new()$gene
    
    # Use character(0) to remove all choices when input hasn't been received
    if (is.null(x))    x <- character(0)
    
    # Update the selection dropdown menu label based on length of input
    if (length(x) == 1)      text_select_input <- " input)"
    else if (length(x) > 1)  text_select_input <- " inputs)"
    else                     text_select_input <- NULL
    
    # Produce dynamically updating dropdown menu to select gene from inputs
    updateSelectInput(session, "pick_timecourse",
                      label = paste("Select gene to display (from ", length(x), text_select_input),
                      choices = x,
                      selected = head(x, 1) # Select first one by default
    )
  })
  
  # --- STATIC TIMECOURSE ---
  
  # Generate ribbon plot and save output for user download
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
    
    # Perform input validation
    validate(
      need(is.null(anno_genes),
           glue("\n\n\nThe input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\nThe input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Check if expression is all zero in the selected brain region
    # TODO: fix this to make it work
    all_zero <- ribbon_plot(gene   = input$pick_timecourse,
                      region = input_new()$region)$zero
    
    validate(
      need(all_zero == FALSE, "This gene has no detected expression in the selected brain region.")
    )
    
    # Create plot with custom plotter (functions.R)
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
  
  # --- INTERACTIVE TIMECOURSE --- 
  
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
    
    # Perform input validation
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
    
    # Create plot with custom plotter (functions.R)
    ribbon_plot(gene   = input$pick_timecourse,
                region = input_new()$region,
                make_plotly = TRUE)$plot
    
  })
  
  # Output the interactive plot
  output$plotlyRibbon <- renderPlotly({ 

      # Position legend to the right of the plot
      layout(ribbon_plotly(), legend = list(x = 1, y = 0))
  })

  # Download static timecourse plot as a PDF
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
    
    # Perform input validation
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
  
  # Create a dimensionality reduction plot coloured by cluster
  # See dr_plot() in functions.R
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
  
  # Store the coordinates for the center (average x & y) of each cluster
  # See dr_plot() in functions.R
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
  
  # Output the hover tooltip UI
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
  
  # Store the gene expression for the gene(s) in user input
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
  # See vln() in functions.R
  output$vln_joint <- renderPlot({
    
    vln(dr_joint_exp(),
        palette = input_new()$clust_palette,
        points  = input_new()$vln_points)
    
  })
  
  #### ---- Sample joint analysis tab content ----
  
  # Store data for the selected region, format dataframe
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
    
    # Map the dr_plot() plotter (functions.R) to each of the timepoints, 
    #   creating one plot per timepoint
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
  
  # Get gene expression for the gene(s) in user input
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
    
    # Map vln() plotter (functions.R) to each timepoint, producing
    #   one plot per timepoint
    map(dr_sample_exp(),
        ~ vln(.x,
              palette = pal,
              points = input_new()$vln_points) +
          theme(plot.margin = unit(c(0.5, 0, 1, 1.5), "cm"))) %>%
      {plot_grid(plotlist = ., ncol = 1, align = "hv",
                 labels = timepoints, label_size = 15)}
    
  })
  
  #### ---- Clusters ranked by expression tab content ----
  
  # Dynamic dropdown menu to choose which input gene to plot
  observe({
    x <- input_new()$gene
    
    # Use character(0) to remove all choices when input hasn't been received
    if (is.null(x))    x <- character(0)
    
    # Update the selection dropdown menu label based on length of input
    if (length(x) == 1)      text_select_input <- " input)"
    else if (length(x) > 1)  text_select_input <- " inputs)"
    else                     text_select_input <- NULL
    
    # Produce dynamically updating dropdown menu to select gene from inputs
    updateSelectInput(session, "pick_ranked",
                      label = paste("Select gene to display (from ", length(x), text_select_input),
                      choices = x,
                      selected = head(x, 1) # Select first one by default
    )
  })
  
  # Store the plot
  ranked_plot <- reactive({
    
    # Check whether a gene was provided or not
    validate(
      need(length(input_new()$gene) > 0, "\n\n\n      Please enter a gene.")
    )
    
    # Check all genes OR selected genes against the dataset & annotations
    #   (need to check all genes if mean expression over genes is being displayed)
    if(input_new()$mean_exp){
      # Check ALL inputs against the dataset & annotations
      num_genes <- NULL
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
    
    # Perform input validation
    validate(
      need(is.null(anno_genes),
           glue("\n\n\n     The input gene \"{anno_genes}\" is in the gene annotation but was not detected in this dataset."))
    )
    validate(
      need(is.null(not_anno_genes),
           glue("\n\n\n     The input gene \"{not_anno_genes}\" was not found in the gene annotation."))
    )
    
    # Customize plot axis labels & title for plotting mean expression vs. one gene
    if (input_new()$mean_exp){
      df <- bubble_prep(gene = input_new()$gene,
                        show_mean = TRUE) %>% 
        filter(Gene == "MEAN")
      y_axis <- "Mean gene expression"
      title <- "Mean expression over all selected genes"
    } else {
      df <- bubble_prep(gene = input$pick_ranked)
      y_axis <- glue("Mean {input$pick_ranked} expression")
      title <- input$pick_ranked
    }
    
    # Create the plot using custom plotter (functions.R)
    ranked_exp_plot(df, title, y_axis)
  })

  # Output the plot
  output$rank_tick_plot <- renderPlot({
    ranked_plot()
  })
  
  # Download plot as a PDF
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
  
  # Store initial dimensions of plot, will show briefly while plot is loading,
  #   so set large height value so the plot that shows briefly is out of view
  #   (Need to initialize this in order to later store dimensions dynamically)
  heatmap_dims <- reactiveValues(height = "100in", width = "12in")
  
  # Store heatmap plot
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
    
    # Perform input validation
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
    
    # Plot heatmap based on gene expression values for input genes,
    # selected cell classes, and selected annotations using
    # custom plotter (see functions.R)
    hm_out <- heatmap_anno_plot(df = bubble_prep(gene = input_new()$gene),
                                cell_classes = input_new()$heatmap_cells,
                                anno = input_new()$heatmap_anno)
    
    # Modify reactive values to dynamically plot the heatmap width
    heatmap_dims$height <- hm_out$height 
    heatmap_dims$width <- hm_out$width
    
    # Output the heatmap
    hm_out$plot
    
  })
  
  output$heatmap <- renderPlot(heatmap_plot())
  
  # Plot the heatmap using dynamic width and height values stored above
  output$heatmapUI <- renderUI({
    plotOutput("heatmap", 
               width = heatmap_dims$width, 
               height = heatmap_dims$height)
    })
  
  # Download heatmap plot as a PDF
  output$download_heatmap <- 
    downloadHandler(filename = "heatmap.pdf",
                    content = function(file) {
                      ggsave(file, 
                             heatmap_plot(),
                             # Adjust pdf dimensions to display the full heatmap
                             width = (as.numeric(substr(heatmap_dims$width,
                                                       1, nchar(heatmap_dims$width)-2)))*1.2,
                             
                             height = (as.numeric(substr(heatmap_dims$height,
                                                        1, nchar(heatmap_dims$height)-2))) + 1.5,
                             units = "in", 
                             scale = 1,
                             limitsize = FALSE # Allow large pdfs
                             )
                    },
                    contentType = "application/pdf")
  
}

