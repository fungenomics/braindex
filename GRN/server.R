#----------------------------------problem---------------------------------
#need to check select input for TF and genes and display any that are not in the data or annotation
#then generate a list of genes or TFs to use and then supply that as the input_new()$tf input

server <- function(input, output, session) {

  # Dynamic UI, change the selectInput tf lists on display depending on the brain region that is selected
  updateSelectizeInput(session, inputId = "TF", choices = all_tf_list, 
                       selected = c("Arx","Lef1"), server = TRUE)
  
  observeEvent(input$region,{
    if(input$region == "cortex"){
       updateSelectizeInput(session, inputId = "TF", choices = all_tf_list, 
                            selected = c("Arx","Lef1"), server = TRUE)
      updateSelectizeInput(session, inputId = "gene", choices = unique(data_cortex_extended$TF_target_gene_info$gene), 
                           selected = c("Dlx6","Sox6"), server = TRUE )
      
    }
    else{
      updateSelectizeInput(session, inputId = "TF", choices = all_tf_list, 
                            selected = c("Lhx5","Pax7"), server = TRUE)
      updateSelectizeInput(session, inputId = "gene", choices = unique(data_pons_extended$TF_target_gene_info$gene), 
                           selected = c("Gad2"), server = TRUE)
      
    }
    #updateRadioButtons(session, "show", selected = "stop") #resets the network visualization 
  })

  
  #uses the input update button to update a list of the parameters of the app for the following functions
  input_new <- eventReactive(input$update,{
    #data_cortex and data_pons are created by data_prep.R and loaded in global.R contains a list 
    #data files (defined in data_prep.R)
    l <- list()
    if(input$region == "cortex"){
      l <- data_cortex_extended
      temp <- paste("F-", input$time, sep = "")
      l$dend_order <- dend_order_joint_cortex_extended
    }
    else if(input$region == "pons"){
      l <- data_pons_extended
      temp <- paste("P-", input$time, sep = "")
      l$dend_order <- dend_order_joint_pons_extended
      }
    #str(input$TF)
    l$tf <- input$TF
    l$region <- input$region
    #l$input_pathway <- input$input_pathway
    l$method <- input$method
    #l$num_cell_plot <- input$num_cell_plot
    l$time_point <- temp
    l$gene <- input$gene
    l$label <- input$label
    l$dim_red <- input$dim_red
    l$cluster_label <- input$cluster_label
    l$grn_toggle <- input$grn_toggle
    l$table_toggle <- input$table_toggle
    l$heatmap_toggle <- input$heatmap_toggle
    l$cluster_toggle <- input$cluster_toggle
    l$as_cluster <- input$active_specific_cluster
    l$as_toggle <- input$as_toggle
    l$fc <- input$fc
    l$as_tp <- input$as_tp
    # l$gene_file_path <- input$file_gene$datapath
    # print(l$gene_file_path)
    # l has following elements with same names for both options above:
    # l contains ...
    
    # We will use the same name attributes to retrieve data
    return (l)
    })
  
  #input_tf <- reactive(input_new()$tf)
  

  #------------------------problem--------------------------
  #once the per sample toggle is hit, there is brief moment where a pathing error 
  #is displayed because the paths depend on the time input for each tab and the time 
  #input is not reactive on the update button, but reaction on use input
  #so the pathing error results from the app trying to access the file using a NULL time value
  #before the insertUI kicks in and supplies the correct time value to the path
  #using priorities in the observeEvent has not helped 
  
  #solved using req() lol
  
  observeEvent(input$tabs, {
    if (input$table_toggle == TRUE & !identical(input$tabs,
                                                "Transcription Factor Target Information")){
      removeUI(
        selector = "div:has(>> #table_tp)"
      )
    }
    else if (input$table_toggle == TRUE & identical(input$tabs,
                                                    "Transcription Factor Target Information")){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "table_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
    }
  }, priority = 1000)
  
  #insert timepoint selection when the toggle for sample data in the table tab is acted on
  observeEvent(input$table_toggle,{
    if(input$table_toggle){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "table_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
      #print(input$tabs)
      output$table_data <- renderText({
        glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$table_tp}</b>")
      })
    }
    else if (input$table_toggle == FALSE){
      removeUI(
        selector = "div:has(>> #table_tp)"
      )
      output$table_data <- renderText({
        ""
      })
    }
  }, priority = 1000) 
  
  observeEvent(input$tabs, {
    if (input$grn_toggle == TRUE & !identical(input$tabs,
                                                "Regulatory Network Visualization")){
      removeUI(
        selector = "div:has(>> #grn_tp)"
      )
    }
    else if (input$grn_toggle == TRUE & identical(input$tabs,
                                                    "Regulatory Network Visualization")){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "grn_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
    }
  })
  #insert timepoint selection when the toggle for sample data in the GRN tab is acted on
  observeEvent(input$grn_toggle,{
    if(input$grn_toggle){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "grn_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
      output$grn_data <- renderText({
        glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$grn_tp}</b>")
      })
    }
    else if (input$grn_toggle == FALSE){
      removeUI(
        selector = "div:has(>> #grn_tp)"
      )
      output$grn_data <- renderText({
        ""
      })
    }
  }) 
  
  observeEvent(input$heatmap_toggle,{
    if(input$heatmap_toggle){
      updateCheckboxGroupInput(session, inputId = "method", "Cluster Method",
                             choices = c("Per Sample Clustering" = "joint"),
                             selected = "joint")
      updateSelectInput(session, inputId = "time",
                                             label = "Time-point to Visualize",
                                             choices = dev_time_points,
                                             selected = "e12")
      output$hm_data <- renderText({
        glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$time}</b>")
      })
    }
    else{
      updateCheckboxGroupInput(session, inputId = "method", "Cluster Method",
                               choices = c("Joint Cluster" = "joint",
                                           "Sample Cluster" = "Cluster"),
                               selected = "joint")
      updateSelectInput(session, inputId = "time",
                        label = "Time-point to Visualize",
                        choices = c("All Time-Points", dev_time_points),
                        selected = "All Time-Points")
      output$hm_data <- renderText({
        ""
      })
    }
    
  })
  
  observeEvent(input$tabs, {
    if (input$cluster_toggle == TRUE & !identical(input$tabs,
                                              "Clustering")){
      removeUI(
        selector = "div:has(>> #cluster_tp)"
      )
    }
    else if (input$cluster_toggle == TRUE & identical(input$tabs,
                                                  "Clustering")){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "cluster_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
    }
  })
  #insert timepoint selection when the toggle for sample data in the GRN tab is acted on
  observeEvent(input$cluster_toggle,{
    if(input$cluster_toggle){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "cluster_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
      output$dr_data <- renderText({
        glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$cluster_tp}</b>")
      })
    }
    else if (input$cluster_toggle == FALSE){
      removeUI(
        selector = "div:has(>> #cluster_tp)"
      )
      output$dr_data <- renderText({
        ""
      })
    }
  })
  
  #Update the active and specific tab inputs and responding to toggle
  observeEvent(input$tabs, {
    if (input$as_toggle == TRUE & !identical(input$tabs,
                                                  "active_specific")){
      removeUI(
        selector = "div:has(>> #as_tp)"
      )
    }
    else if (input$as_toggle == TRUE & identical(input$tabs,
                                                      "active_specific")){
      #print("hello")
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "as_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
    }
  })
  #insert timepoint selection when the toggle for sample data in the GRN tab is acted on
  observeEvent(input$as_toggle,{
    if(input$as_toggle){
      insertUI(
        selector = "#region",
        where = "afterEnd",
        ui = selectInput(inputId = "as_tp",
                         label = "Developmental Time-Point to Explore",
                         choices = dev_time_points,
                         multiple = FALSE,
                         selected = "e12")
      )
      output$as_data <- renderText({
        glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$as_tp}</b>")
      })
    }
    else if (input$as_toggle == FALSE){
      removeUI(
        selector = "div:has(>> #as_tp)"
      )
      output$as_data <- renderText({
        ""
      })
    }
  })
  
  # observeEvent(input$tabs, {
  #   if (input$bb_toggle == TRUE & !identical(input$tabs,
  #                                             "bubble")){
  #     removeUI(
  #       selector = "div:has(>> #bb_tp)"
  #     )
  #   }
  #   else if (input$bb_toggle == TRUE & identical(input$tabs,
  #                                                 "bubble")){
  #     insertUI(
  #       selector = "#region",
  #       where = "afterEnd",
  #       ui = selectInput(inputId = "bb_tp",
  #                        label = "Developmental Time-Point to Explore",
  #                        choices = dev_time_points,
  #                        multiple = FALSE,
  #                        selected = "e12")
  #     )
  #   }
  # })
  # #insert timepoint selection when the toggle for sample data in the GRN tab is acted on
  # observeEvent(input$bb_toggle,{
  #   if(input$bb_toggle){
  #     insertUI(
  #       selector = "#region",
  #       where = "afterEnd",
  #       ui = selectInput(inputId = "bb_tp",
  #                        label = "Developmental Time-Point to Explore",
  #                        choices = dev_time_points,
  #                        multiple = FALSE,
  #                        selected = "e12")
  #     )
  #     output$bb_data <- renderText({
  #       glue("<b>Current Dataset: {str_to_title(input_new()$region)} - Time-Point {input$bb_tp}</b>")
  #     })
  #   }
  #   else if (input$bb_toggle == FALSE){
  #     removeUI(
  #       selector = "div:has(>> #bb_tp)"
  #     )
  #     output$bb_data <- renderText({
  #       ""
  #     })
  #   }
  # }) 
  
  #update inputs when toggle is turned so that the plots auto-update
  observeEvent(input$grn_toggle|input$table_toggle|input$heatmap_toggle|input$cluster_toggle|input$as_toggle|input$bb_toggle, {
    click(id = "update")
  }, ignoreInit = TRUE, priority = 1)
  
  #updates a reactive value reg depending on the input region which is used to select the right dataset to display in the app
  #uses the input_new() region because wants to be dependant on the update button
  reg <- reactive({
    if(identical(input_new()$region, "cortex")){
      "ct"
    }
    else{
      "po"
    }
  })
  
  observeEvent(input$info|input$info1|input$info2|input$info3|input$info4, {
    sendSweetAlert(session, title = "What Is This?", text = 
                     "Use toggle to explore data from SCENIC analyses performed
                   on each developmental time-point individually.")
  }, ignoreInit = TRUE)
  
  tf_list <- reactiveValues(
    TF_in_data = NULL,
    TF_not_data = NULL
  )
  output$tf_check <- renderText({
    # if(identical(input$tabs, "active_specific")){
    #   ""
    # }
    if(!length(tf_list$TF_in_data) > 0 & !is.null(tf_list$TF_in_data)){
      "<font color=\"#FF0000\"><b> None of the input TFs are active in the current dataset. </b></font>"
    }
    else{
      glue("The following TFs in your input are not active in the dataset
         that your are currently exploring: <font color=\"#FF0000\"><b> {toString(tf_list$TF_not_data)} </b></font>")
    }
  })
  
  # -----------------------------Table------------------------------------------

  #filter the data, add a column for logos, then display
    output$table1 <- renderDataTable({
        # process data, filter the lines with our interested TF
        
      datafile <- glue("data_{reg()}_{input$table_tp}")
      
      
      if(input$table_toggle){
        
        req(input$table_tp)
        
        temp <- check_tf_input(input_new()$tf, get(datafile)$unique_TF)
        
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        #print(tf_list)
        
        subset_data <- get(datafile)$TF_target_gene_info %>% 
          dplyr::filter(TF %in% tf_list$TF_in_data) %>% 
          select(TF, gene, starts_with("Genie3Weight"), highConfAnnot, nMotifs, bestMotif)
        
        #the column name for the SCENIC inference weight has a different name depending on if the timepoint
        #is in the original dataset from 2019 nat genetic or if it is the extended dataset
        #so need to have a variable weight_col to rename the appropriate column in the datatable() function below
        if(input$table_tp %in% c("e12", "e15", "p0", "p3", "p6")){weight_col <- "Genie3Weight.weight"}
        else{weight_col <- "Genie3Weight"}
           
      }
      
      else{
        temp <- check_tf_input(input_new()$tf, input_new()$unique_active_TFs_bare)
        
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data

        #print(tf_list)
        
        subset_data <- input_new()$TF_target_gene_info %>% 
          dplyr::filter(TF %in% tf_list$TF_in_data) %>% 
          select(TF, gene, starts_with("Genie3Weight"), highConfAnnot, nMotifs, bestMotif)
        weight_col <- "Genie3Weight"
      }
      
      subset_data <- addMotifPic(subset_data)
      
      datatable(subset_data, escape = FALSE,
                colnames = c('Gene' = 'gene', 'Number of Motifs' = 'nMotifs',
                             'Best Motif' = 'bestMotif', 
                             'Strength of Association' = weight_col,
                             'Logo' = 'motif_logo'),
                rownames = FALSE)
    })
  

  # -----------------------------GRN------------------------------------------
    gene_list <- reactiveValues( #makes reactive values to take in the user input gene list
      data = NULL,
      clear = FALSE
    )
    
    observe({ #reads in the gene_list into a vector form if the clear flag is false
      req(input$file_gene)
      req(!gene_list$clear)

      gene_list$data <- read_csv(input$file_gene$datapath) %>% deframe() %>% map_chr(str_to_title) 
      #reads the input file and turns it into a vector and then capitalizes each gene name
      
    })
    
    observeEvent(input$file_gene, { #whenever a file is inputted, independent of the update click, clear is set to 
      #false so the data is read in
      gene_list$clear <- FALSE
    }, priority = 1000)
    
    observeEvent(input$reset, { #once the reset button is pressed, the data is removed and clear is set to true 
      #so file is not read in, all of this independent of the update master button that is usually used 
      gene_list$data <- NULL
      gene_list$clear <- TRUE
    }, priority = 1000)
    
    output$file1_ui <- renderUI({
      input$reset ## Create a dependency with the reset button
      fileInput(
        "file_gene", "Choose a CSV file containing your genes list",
        accept = c(
          "text/csv",
          "text/comma-separated-values, text/plain",
          ".csv"),
        multiple = FALSE,
        placeholder = "example_list.csv"
      )
    })
    #check if there is a user input gene_list file, if there is, use it, if not, use the selectInput genes 
    igraph_network <- reactive ({
      
      if(is.null(gene_list$data)){
        gene_to_highlight <- input_new()$gene 
      }
      
      else{
        gene_to_highlight <- gene_list$data
      }
      
      datafile <- glue("data_{reg()}_{input$grn_tp}")
      
      if(input$grn_toggle){
        req(input$grn_tp)
        
        temp <- check_tf_input(input_new()$tf, get(datafile)$unique_TF)
        
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        req(length(tf_list$TF_in_data) > 0)
        
        make_network(tf_list$TF_in_data, get(datafile)$TF_target_gene_info,
                     gene_to_highlight)
      }
      
      else{
        
        temp <- check_tf_input(input_new()$tf, input_new()$unique_active_TFs_bare)
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        req(length(tf_list$TF_in_data) > 0)
         
        make_network(tf_list$TF_in_data, input_new()$TF_target_gene_info,
                     gene_to_highlight) #returns an igraph network object
      }
    })
    network_ggplot <- reactive({
      plot_network(igraph_network(), input_new()$label, tf_list$TF_in_data)
    })
    
    output$network <- renderPlotly({
      net_plotly <- network_ggplot() %>% ggplotly(height = 700, tooltip = "text") %>% 
        layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
               hovermode = "x", hoverdistance = 100)
      net_plotly
      
    })
    
    output$download_network <- downloadHandler(filename = "network.pdf",
                                                contentType = "application/pdf",
                                                content = function(file){
                                                  ggsave(filename = file, plot = network_ggplot(),
                                                         width = 8.5, height = 11)
                                                })
    
    
    # -----------------------------Heatmap-------------------------------------------
    # output$color_hm_palette <- renderImage({
    #   
    #   expr = list(src = "www/timeseries_color.png", #picture of colors corresponding with clusters 
    #        alt = "This is alternate text")
    #   
    #   
    # },
    # deleteFile = FALSE)
    
    hm_joint_cluster_plot <- reactive({
      
      req("joint" %in% input_new()$method)
      
      if(input_new()$heatmap_toggle == TRUE){
        
        datafile <- glue("data_{reg()}_{input$time}")
        
        temp <- check_tf_input(input_new()$tf, unique(get(datafile)$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
         
        req(length(tf_list$TF_in_data) > 0)
        
        req(input$time != "All")
        
        plot_heatmap(tf_list$TF_in_data, "joint",input_new()$region,
                     get(datafile)$TF_and_ext, per_sample = TRUE,
                     timepoint = input$time)
      }
      else{
        temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        req(length(tf_list$TF_in_data) > 0)
        
        plot_heatmap(tf_list$TF_in_data, "joint",input_new()$region,
                   input_new()$TF_and_ext)
      }

    })
    #need to adjust this, at least zoomable, change color scheme to dark blue and neon yellow
    hm_sample_cluster_plot <- reactive({ #this is displaying heatmap clustering by sample cluster and not joint cluster
      req("Cluster" %in% input_new()$method)
      
      temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
      tf_list$TF_in_data <- temp$TF_in_data
      tf_list$TF_not_data <- temp$TF_not_data
      
      req(length(tf_list$TF_in_data) > 0)
      
      plot_heatmap(tf_list$TF_in_data, "Cluster",input_new()$region, 
                   input_new()$TF_and_ext, timepoint = input_new()$time_point)
      
    })
 
    output$heatmap_joint <- renderPlot({
      #print(str(input_new()$tf))
      req(length(input_new()$tf) != 23 ) #change this line after the tf input check has been done
      hm_joint_cluster_plot()
    })
    
    hm_name <- reactive({
      if(input$heatmap_toggle){
        glue("{reg()}_{input$time}_heatmap.pdf")
      }
      else{
        "heatmap_joint.pdf"
      }
    })
    output$download_hm_joint <- downloadHandler(filename = hm_name(),
                                               contentType = "application/pdf",
                                               content = function(file){
                                                 ggsave(filename = file, plot = hm_joint_cluster_plot(),
                                                        width = 20, height = 25)
                                               })
    
    output$heatmap_cluster <- renderUI({
      
      if(input$time == "All"){plot_width = '3500px'}
      else{plot_width = '800px'}
      
      plotOutput('heatmap_cluster_plot', width = plot_width)             
    })
    
     output$heatmap_cluster_plot <- renderPlot({
       
      hm_sample_cluster_plot() 
      

    })
    
    output$download_hm_cluster <- downloadHandler(filename = "heatmap_cluster.pdf",
                                               contentType = "application/pdf",
                                               content = function(file){
                                                 ggsave(filename = file, plot = hm_sample_cluster_plot(),
                                                        width = 20, height = 25)
                                               })
 
    
    # -----------------------------DR plots------------------------------------------ 
    # The cluster scatterplot is always plot by cells, so we use an independent reactive
    # value for this plot
    activity_data_cluster <- reactive({
      
      
      if(input_new()$cluster_toggle){
        
        datafile <- glue("data_{reg()}_{input$cluster_tp}")
        #print(datafile)
        
        temp <- check_tf_input(input_new()$tf, unique(get(datafile)$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        req(length(tf_list$TF_in_data) > 0)
        
        req(input$cluster_tp)
        
        create_activity_data(tf_list$TF_in_data, "Cell", input_new()$region,
                             get(datafile)$TF_and_ext, per_sample = input_new()$cluster_toggle,
                             timepoint = input$cluster_tp, bad_cells = get(datafile)$bad_cells)
      }
      else{
        temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data
        tf_list$TF_not_data <- temp$TF_not_data
        
        req(length(tf_list$TF_in_data) > 0)
        # use the feature of feather data to read certain col to optimize speed
        create_activity_data(tf_list$TF_in_data, "Cell",input_new()$region, input_new()$TF_and_ext)
      }
    }) #this part takes the TF that the user selects and subsets the feather file containing all 
    #regulon activity data to include only the regulons that the user selects and this is entered into the 
    #UMAP plot function 
    


    Umap_plot_1 <- reactive({
      req(length(input_new()$tf)>0)
      
      if(input_new()$cluster_toggle){
        
        datafile <- glue("data_{reg()}_{input$cluster_tp}")
        
        req(input$cluster_tp)
        
        plot_dr(tf_number = 1, get(datafile)$cell_metadata, 
                  activity_data_cluster(), input_new()$dim_red)
      }
      else{
        plot_dr(tf_number = 1,input_new()$cell_metadata, 
                  activity_data_cluster(), input_new()$dim_red)
      }
    })
    Umap_plot_2 <- reactive({
      req(ncol(activity_data_cluster()) > 2)
      
      if(input_new()$cluster_toggle){
        
        datafile <- glue("data_{reg()}_{input$cluster_tp}")
        
        req(input$cluster_tp)
        
        plot_dr(tf_number = 2, get(datafile)$cell_metadata, 
                  activity_data_cluster(), input_new()$dim_red)
      }
      else{
        plot_dr(tf_number = 2,input_new()$cell_metadata,
                  activity_data_cluster(), input_new()$dim_red)
      }
    })
    
    output$color_by_cluster <- renderPlot({
      if (input_new()$cluster_toggle){
        
        datafile <- glue("data_{reg()}_{input$cluster_tp}")
        
        req(input$cluster_tp)
        
        color_by_cluster(get(datafile)$cell_metadata, 
                         input_new()$dim_red, input_new()$cluster_label,
                         per_sample = input_new()$cluster_toggle)
      }
      else{
        color_by_cluster(input_new()$cell_metadata, 
                         input_new()$dim_red, input_new()$cluster_label)
      }
    })
    
    output$cluster1 <- renderPlot({
      Umap_plot_1()
    })
    output$download_UMAP_1 <- downloadHandler(filename = "UMAP1.pdf",
                                              contentType = "application/pdf",
                                              content = function(file){
                                                ggsave(filename = file, plot = Umap_plot_1(),
                                                       width = 20, height = 20)
                                              })
    output$cluster2 <- renderPlot({
      Umap_plot_2()

    })
    
   
    output$download_UMAP_2 <- downloadHandler(filename = "UMAP2.pdf",
                                              contentType = "application/pdf",
                                              content = function(file){
                                                ggsave(filename = file, plot = Umap_plot_2(),
                                                       width = 20, height = 20)
                                              })
  
    
    
    
    # --------------------------------------Timeseries-------------------------------------------
    # tf_nexist_data <- reactive({
    #   tf_nexist <- ""
    #   for(tf in input_new()$tf){
    #     if (tf %in% input_new()$tfs_not_exist_timeseries){
    #       tf_nexist <- paste(tf_nexist,tf,sep = " ")
    #     }
    #   }
    # })
    # tf_desc_timeseries <- reactive({
    #   tf_nexist_string <- ""
    #   for(tf_n in input_new()$tfs_not_exist_timeseries){
    #     tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
    #   }
    #   text <- glue("We do not have these followning tfs in this tab: {tf_nexist_string}")
    #   
    #   tf_nexist <- ""
    #   for(tf in input_new()$tf){
    #     if (tf %in% input_new()$tfs_not_exist_timeseries){
    #       tf_nexist <- paste(tf_nexist,tf,sep = " ")
    #     }
    #   }
    # 
    #   if(tf_nexist == ""){
    #     text <- "Good! All of your input tfs exist in our timeseries activity datasets!"
    #   }
    #   else{
    #     tf_nexist_string <- ""
    #     for(tf_n in input_new()$tfs_not_exist_timeseries){
    #       tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
    #     }
    #     text <- glue('Those tfs in your input list does not not exist in our
    #                timeseries datasets: {tf_nexist}.
    #                We do not have these followning tfs in this tab: {tf_nexist_string}')
    #   }
    # })
    
    # output$tf_timeseries_desc <- renderText({
    #   # tf_desc_timeseries()
    #   tf_nexist_string <- ""
    #   for(tf_n in input_new()$tfs_not_exist_timeseries){
    #     tf_nexist_string <- paste(tf_nexist_string,tf_n,sep = " " )
    #   }
    #   text <- glue("We do not have data for the following trancription: {tf_nexist_string}")
    #   
    #   
    # })
    
    
    output$timeseries_desc <- renderText({
      text <- "Click option: double clicking a cell type in the legend displays that cell type ONLY;
      single click removes that cell type from the plot. Mouse over ribbons in the plot to see the cell types. 
      We only support four plots of your first four transcripton factor inputs."
    
    })
    
    # we must transform the TF format, from raw form (Arx) to (Arx_extended (21g)) to fetch
    # information
    TF_transformed <- reactive({
      
      temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
      tf_list$TF_in_data <- temp$TF_in_data
      tf_list$TF_not_data <- temp$TF_not_data
      
      req(length(tf_list$TF_in_data) > 0)
      
      translate_tf(input_new()$tf,input_new()$binary_active_TFs)
      })
    
    # ggplotly_list_plot <- reactive({
    #   req(TF_transformed())
    #   # binary_active_TFs is loaded at beginning by data_prep.R
    #   plot_list <- lapply(TF_transformed(), plot_timeseries, cell_metadata = data_cortex$timeseries_input_meta, 
    #                       activity = data_cortex$binary_activity, make_plotly = TRUE)
    #   # produce a list of ggplotly plots
    #   subplot(plot_list, nrows = 2, margin = 0.04, heights = c(0.6, 0.4), shareX = TRUE, shareY = FALSE)
    #   
    # })
    
    
    ggplot_list_plot <- reactive({
      req(TF_transformed())
      plot_list <- lapply(TF_transformed(), plot_timeseries, cell_metadata = input_new()$timeseries_input_meta, 
                          activity = input_new()$binary_activity, make_plotly = FALSE, show_legend = FALSE)
      plot_grid(plotlist = plot_list)
    })
    
    output$timeseries1 <- renderPlotly({ # a plotly list
      req(length(input_new()$tf)>0)
      plot_timeseries(TF_transformed()[1][1], input_new()$timeseries_input_meta, input_new()$binary_activity, make_plotly = TRUE)
     })
    
    output$timeseries2 <- renderPlot({ # a ggplot list
      ggplot_list_plot()
      
    })
    
    output$download_ribbon_1 <- downloadHandler(filename = "timeseries_ribbon.pdf",
                                                contentType = "application/pdf",
                                                content = function(file){
                                                  ggsave(filename = file, plot = ggplot_list_plot(),
                                                         width = 20, height = 15)
                                                })
    
    # output$pons_timeseries <- renderImage({
    #      
    #      list(src = "www/pons_timeseries.png",
    #           alt = "This is alternate text")
    #      
    #    }, deleteFile = FALSE)
    
    output$cell_proportion_timeseries <- renderPlotly({
      
      if(identical(input_new()$region, "cortex")){
        ggplotly(timeseries_proportion_plots$forebrain_plot + 
                   ggtitle("Proportion of cells from each major cell class in the forebrain over the time course."), tooltip = "Cluster") %>% 
          style(hoveron = "points + fills")
      }
      else{
        ggplotly(timeseries_proportion_plots$pons_plot + 
                   ggtitle("Proportion of cells from each major cell class in the pons over the time course."), tooltip = "Cluster") %>% 
          style(hoveron = "points + fills")
      }
      
    })
   
    
    # output$timeseries_color <- renderImage({
    #   
    #   list(src = "www/timeseries_color.png",
    #        alt = "This is alternate text")
    #   
    # },
    # deleteFile = FALSE)
    
#-------------------------------------Active specific-------------------
    #same idea as reg but used to update the cluster list in this tab, uses the input$region because do not want 
    #dependence on the update button
    reg2 <- reactive({
      if(identical(input$region, "cortex")){
        "ct"
      }
      else{
        "po"
      }
    })
    clust_list <- reactive({
      if(input$as_toggle){
        
        req(input$as_tp)
        
        datafile <- glue("data_{reg2()}_{input$as_tp}")
        get(datafile)$cell_metadata %>% select(ID_20201028_with_exclude) %>% unique() %>%
          filter(!grepl("EXCLUDE", ID_20201028_with_exclude)) %>% deframe()
        
      }
      else{
        datafile <- glue("data_{input$region}_extended")
        
        get(datafile)$cell_metadata %>% select(Sample_cluster) %>% unique() %>% 
          filter(!grepl("EXCLUDE", Sample_cluster)) %>% deframe()
      }
    })
    update_in <- observe({
      updateSelectizeInput(session, inputId = "active_specific_cluster", choices = clust_list(), 
                           selected = clust_list()[1], server = TRUE)
    }, priority = 1000)
    
    active_specific_data <- reactive({
      if(input$as_toggle){
        
        req(input$as_tp)
        req(grepl(input$as_tp, input_new()$as_cluster))
        
        data_sample <- glue("{reg()}_{input$as_tp}")
        sample_name <- glue("data_{reg()}_{input$as_tp}")
      }
      else{
        data_sample <- glue("joint_{input_new()$region}_extended")
        sample_name <- glue("data_{input_new()$region}_extended")
      }

      #print(input_new()$as_cluster)
      woof <- list(
        "data" = active_specific_prep(data_sample, input_new()$as_cluster),
        "name" = sample_name
      )
      #active_specific_prep(data_sample, input_new()$as_cluster)
    })
    
    
    output$as_clust <- renderUI({
       fluidRow(
         column(width = 7, plotOutput("active_specific_scatter")),
         column(width = 5, 
                
                (div(style='height:400px;overflow-y: scroll; width:450px;overflow-x: scroll;',
                     tableOutput("active_specific_table"))))
       )
      #fluidRow(plotOutput("active_specific_cluster"))
    })
    
    output$as_tf <- renderUI({
      plotOutput("as_bar_AUC", height = '800px')
    })
    
    output$active_specific_dr <- renderPlot({
      data <- get(active_specific_data()$name)
      #cluster <- gsub(".+_", "", input_new()$as_cluster)
      awo <- hm_anno$side_colors$Cluster
      names(awo) <- gsub(" ", "_", names(hm_anno$side_colors$Cluster))
      to_color <- replace(awo, !grepl(input_new()$as_cluster, names(awo)), "#e5e5e5")
      #print(to_color)
      gg <- color_by_cluster(data$cell_metadata, 
                             "umap", FALSE,
                             per_sample = input_new()$as_toggle)
      if(input$as_toggle == FALSE){
        gg <- gg + geom_point(aes(color = Sample_cluster), alpha = 0.2)
      }
      gg <- gg + scale_color_manual(values = to_color) + 
        theme(legend.position = "none", plot.title = element_text(size = 14, face="bold")) + 
        ggtitle(input_new()$as_cluster)
      gg
    })
    
    #enveloping the bar plot as an reactive expression so it can be called in the download handler 
    as_bar <- reactive({
      if(input$as_toggle == TRUE){
        
        datafile <- glue("data_{reg()}_{input$as_tp}")
        
        req(input$as_tp)
        
        temp <- check_tf_input(input_new()$tf, unique(get(datafile)$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data %>% transform_tf_input(get(datafile)$TF_and_ext) %>% 
          head(4)
        tf_list$TF_not_data <- temp$TF_not_data
        
    
      }
      else{
        temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
        tf_list$TF_in_data <- temp$TF_in_data %>% transform_tf_input(input_new()$TF_and_ext) %>%
          head(4)
        tf_list$TF_not_data <- temp$TF_not_data
        
      }

      req(length(tf_list$TF_in_data) > 0)
      plot_bar_list(active_specific_data()$data$AUC_df, tf_list$TF_in_data)
    }) 
    
    output$as_bar_AUC <- renderPlot({
      as_bar()
    })
    
    
    output$as_bar_download <- downloadHandler(filename = "bar_plot.pdf",
                                                contentType = "application/pdf",
                                                content = function(file){
                                                  ggsave(filename = file, plot = as_bar())                                                })
    
    output$active_specific_scatter <- renderPlot({
      plot_scatter(active_specific_data()$data$tf_table, input_new()$fc, input_new()$as_cluster)
    })
    
    
    output$as_scatter_download <- downloadHandler(filename = "scatter_plot.pdf",
                                              contentType = "application/pdf",
                                              content = function(file){
                                                ggsave(filename = file, plot = plot_scatter(active_specific_data()$data$tf_table, input_new()$fc, input_new()$as_cluster))#why doesnt this work :(
                                              })
    
    output$active_specific_table <- renderTable({
      
      data <- active_specific_data()$data$tf_table %>% select(-is_ext) %>%
        filter(AUC_FC > input_new()$fc) %>% arrange(desc(AUC_in))
      
      temp_col <- data %>% select(TF)
      
      data <- data %>% select(-TF) %>% 
        round(3) %>% mutate(TF = temp_col) %>% 
        transmute('TF' = TF,
                  'Average AUC in Cluster' = AUC_in,
                  'Average AUC in Other' = AUC_out,
                  'AUC Fold Change' = AUC_FC )
      
      # datatable(data, escape = TRUE,
      #           colnames = c('Average AUC in Cluster' = 'AUC_in', 
      #                        'Average AUC in Other' = 'AUC_out',
      #                        'AUC Fold Change' = 'AUC_FC'),
      #           rownames = FALSE)
    })
    
    
    
#-------------------------Bubbles--------------------
    output$dend_image <- renderUI({
      image_source <- switch(input$region,
                             "cortex" = "joint_cortex_extended_tree.png",
                             "pons" = "joint_pons_extended_tree.png")
      image_height <- switch(input$region,
                             "cortex" = "143",
                             "pons" = "120")
      div(style = "margin-top: 3em; margin-bottom: -2em !important;",
          fluidRow(tags$img(src = image_source, width = "1150", height = image_height))
      )
    })
    
    
    #reactively changes data for bubble plot depending on inputs
    bubble_data <- reactive({
      # if(identical(input_new()$region, "cortex")){
      #   dend_source <- dend_order_forebrain_tp
      # }
      # else{
      #   dend_source <- dend_order_pons_tp
      # }
        
      temp <- check_tf_input(input_new()$tf, unique(input_new()$TF_and_ext[["type"]]))
      tf_list$TF_in_data <- temp$TF_in_data %>% transform_tf_input(input_new()$TF_and_ext) %>%
        head(20)
      tf_list$TF_not_data <- temp$TF_not_data
        
      data_sample <- glue("joint_{input_new()$region}_extended")

      #dend_order <- dend_source[data_sample]
      
      bubble_prep(sample = data_sample,
                  tf = tf_list$TF_in_data,
                  dend_order = input_new()$dend_order,
                  scale = FALSE)
      #active_specific_prep(data_sample, input_new()$as_cluster)
    })
    
    num_of_tf <- reactive({
       length(tf_list$TF_in_data)
    })
    
    # Generate the bubbleplot
    output$bubble <- renderPlot({
      
      #print(bubble_data()$data)
      
      plot_bubble(data = bubble_data()$data,
                  label_palette = bubble_data()$label_palette, 
                  input_new()$dend_order)$plot # Get plot part of output
     
    },
    
    # Choose width to align horizontally with dendrogram image
    width = 1143,
    
    # Customize the height of the bubbleplot to scale with the number of genes which
    # are being displayed, after allocating a baseline height for the x-axis & legend
    
    
    height = function() 150 + 30 * num_of_tf()
    )
    
    
    # Create a tooltip with cluster / expression information 
    # that appears when hovering over a bubble 
    # This was adapted from this example: https://gitlab.com/snippets/16220
   
     output$bubble_hover_info <- renderUI({
      
      hover <- input$bubble_hover
      
      
      
      # Find the nearest data point to the mouse hover position
      point <- nearPoints(bubble_data()$data,
                          hover,
                          xvar = "Cluster",
                          yvar = "TF_padded",
                          maxpoints = 1) %>% 
        select(TF, Cluster, AUC, FC)
      
      # Hide the tooltip if mouse is not hovering over a bubble
      if (nrow(point) == 0) return(NULL)
      
      # Create style property for tooltip
      # background is set to the cluster colour, with opacity = 95% ("F2" at end of hex)
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: ", bubble_data()$label_palette[point$Cluster], "F2;",
                      "left: -350px; top: 400px; width: 350px;")
      
      # Set text to white if the background colour is dark, else it's black (default)
      if (dark(bubble_data()$label_palette[point$Cluster])) {
        style <- paste0(style, "color: #FFFFFF")
      }
      
      # Specify text content of tooltips - special content for mean expression plot
        tooltip_text <- paste0("<b> TF: </b>",       point$TF, "<br/>",
                               "<b> Cluster: </b>",    point$Cluster, "<br/>",
                               "<b> TF Activity: </b>",  point$AUC %>% round(3), "<br/>",
                               "<b> TF Activity Fold Change: </b>",     point$FC %>% round(3), "<br/>")
      
      # Actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(tooltip_text))
      )
    })
    
    # Render the bubble plot gene labels separately with ggdraw
    output$bubble_labels <- renderPlot({
      
      ggdraw(plot_bubble(data = bubble_data()$data,
                         label_palette = bubble_data()$label_palette)$labels) # Get labels part of output
      
    },
    
    # Set height of bubble plot gene labels to (hopefully) align with plots
    height = function() 8 + 29 * num_of_tf(),
    
    # Max length of a gene is 200px
    # NOTE: If altering this, also change the corresponding cellWidth for 
    # splitLayout in ui.R
    width = 200
    
    )
    
}

