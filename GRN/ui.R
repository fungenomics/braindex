ui <- fluidPage(
  introjsUI(),
  #useShinyjs(),
  # Application title
  introBox(
    titlePanel("Joint Cortex and Pons Transcription Factor Activity"),
    data.step = 1,
    data.intro = "This app displays transcription factor activity inference data from 
    a developmental timecourse of the mouse Pons and Forebrain."
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # choose which datasets to analyze for the whole app
      actionButton("help", label = "See Instructions"),
      radioButtons("region", "Brain region",
                   # use the names in the vector to display
                   # use the character "joint_cortex" to match the path to import data
                   choices = c("Forebrain" = "cortex",
                               "Pons" = "pons"),
                   selected = "cortex"),
      
      #actionButton("update_tf", label = "Update transcription factors to see the plots"),
      selectInput(inputId = "TF",
                  label = "Transcription Factor",
                  choices = data_cortex$unique_active_TFs_bare,
                  multiple = TRUE,
                  selected = c("Arx","Lef1")),
      # fileInput("file_tf", "Choose CSV File containing your tf list",
      #           accept = c(
      #             "text/csv",
      #             "text/comma-separated-values,text/plain",
      #             ".csv")
      # ),
      
      # 1. table and network graph of related TF and genes
      conditionalPanel(
          condition = "input.tabs == 'Transcription Factor Target Information'",
          
      ),
      conditionalPanel(
          condition = "input.tabs == 'Regulatory Network Visualization'",
          selectInput(inputId = "gene",
                      label = "Genes of Interest",
                      choices = unique(data_cortex$TF_target_gene_info$gene),
                      multiple = TRUE,
                      selected = c("Dlx6","Sox6")),
          #allows user to input a file containing a list of genes to query
          #this input returns a dataframe with 4 columns: name, size, type, datapath
          #each file entered is one row
          
          uiOutput('file1_ui'), ## instead of fileInput('file1', label = NULL) so that the file can be reset 
          # fileInput(
          #   "file_gene", "Choose a CSV file containing your genes list",
          #   accept = c(
          #     "text/csv",
          #     "text/comma-separated-values, text/plain",
          #     ".csv"),
          #   multiple = FALSE,
          #   placeholder = "example_list.csv"
          # ),
          actionButton("reset", label = "Reset File"),
                      
          checkboxInput(inputId = "label", label = "Label Gene Target Nodes", value = FALSE)
         
      ),
      # conditionalPanel(condition = "input.tabs == 'Table and Network'",
      #                 radioButtons("show", "Node Display Option",
      #                               # use the names in the vector to display
      #                               # use the character "joint_cortex" to match the path to import data
      #                               choices = c("Show All Nodes" = "all",
      #                                           #"color by user's input tf list" = "pathway",
      #                                           "Shrink Grey Nodes" = "shrink",
      #                                           "Neglect Grey Nodes" = "neglect",
      # 
      #                                           "Show No Nodes" = "stop"),
      #                               selected = "stop"),#,
      #                  checkboxInput("show_pathway","Color by Input Genes",
      #                                TRUE),
      #                  selectInput("input_pathway", "Gene Pathway of Interest",
      #                              choices = data_cortex$unique_active_TFs_bare,
      #                              multiple = TRUE,
      #                              selected = c("Arx","Lef1"))
      #                  # fileInput("file_gene", "Choose CSV File containing your genes list",
      #                  #           accept = c(
      #                  #             "text/csv",
      #                  #             "text/comma-separated-values,text/plain",
      #                  #             ".csv")
      #                  # )
      #                  #actionButton("update_graph", label = "See the network graph")
      # 
      #  ),
      # 2. heatmap and clustering
      conditionalPanel(condition = "input.tabs == 'Heatmap'",
                       #numericInput(inputId = "num_cell_plot", label = "Number of Cells to Visualize",
                                    #value = 300),
                       # numericInput(inputId = "num_cluster_plot", label = "number of clusters to visualize",
                       #              value = 50),
                       checkboxGroupInput("method", "Cluster Method",
                                          choices = c("Joint Cluster" = "joint",
                                                      "Sample Cluster" = "Cluster"),
                                          selected = "joint"),
                       selectInput(inputId = "time",
                                   label = "Timepoint to Visualize",
                                   choices = c("All","e12", "e15", "p0", "p3", "p6"),
                                   multiple = FALSE,
                                   selected = "e12")
                                          
                       ),
      # 3. time series plot
      conditionalPanel(condition = "input.tabs == 'Time Series'"),
      
      introBox(
      # Update everything
      actionButton("update", label = "Update"),
      data.hint = "click me to update everything!",
      data.step = 2,
      data.intro = "click it to update everything! Do this after you changed your 
      transcription factor input and options. Feel free to QUIT the intro first and update it to see the 
      table and plots",
      data.position = "right"
      ),
    ),
    mainPanel(
      tabsetPanel(
      
        tabPanel(
          p("This tab displays information corresponding to the selected transcription factors and their predicted gene targets."),
          p("- Strength of Association represents "),
          p("- The top row is colored by cluster and the bottom row is colored by transcription factor activity"),
          p("- Only the first two transcriptions factors are displayed"),
          title = "Transcription Factor Target Information",
          textOutput("general_desc"),
          introBox(
            dataTableOutput("table"),
            data.step = 3,
            data.intro = "This table displays the gene targets of the selected transcription factors
              along with information about the genes."
          ),
          value = "Transcription Factor Target Information"
        ),
        
        tabPanel(
          title = "Regulatory Network Visualization",
          #textOutput("general_desc"), # this line breaks things/ probably cause you can't have 2 general_desc
          textOutput("desc"),
          plotOutput("network"),
          #need to work on visualization with the ggNet package
          introBox(
            data.step = 4,
            data.intro = "This table displays a network of your selected transcription factors and
              their top gene targets."
          ),
          value = "Regulatory Network Visualization"
        ),
        # tabPanel(
        #   title = "Table and Network",
        #   textOutput("general_desc"),
        #   introBox(
        # 
        #   dataTableOutput("table"),
        #   data.step = 3,
        #   data.intro = "Table and network tab:
        # A table of tf and its target gene with motifs and other information"
        #   ),
        #   introBox(
        #     data.intro = "Feel free to quit the intro now, click the 'show all nodes' button
        #     in the sidebar to see the cytoscape network graph, then we continue",
        #     data.step = 4
        #   ),
        # 
        # 
        #   textOutput("desc"),
        #   tags$style(type="text/css", "#desc {white-space: pre-wrap;}"),
        #   introBox(
        #   rcytoscapejsOutput("network", width = "1200px",height = "600px"),
        #   data.step = 5,
        #   data.intro = "a network graph visualization displaying detailed information with node color
        #   and size:
        #   Orange nodes are active transcription factors (tf genes that express their own tf);
        #   Purple nodes in the center are your input transcription factors;
        #   Green nodes are your input genes related to input tfs(purple nodes)
        #   ;  grey nodes are other genes."
        #   ),
        #   value = "Table and Network"
        # ),
        
        
        
      tabPanel(
        title = "Heatmap",
        value = "Heatmap",
        fluidRow(
          plotOutput("heatmap_joint")
        ),
        downloadButton("download_hm_joint", "Heatmap by Joint Cluster (PNG)"),
        div(style = "margin-left: 1.3em; margin-right: 1.3em;",
        fluidRow(
          plotOutput("heatmap_cluster")
        )),
        downloadButton("download_hm_cluster", "Heatmap by Timepoint Cluster (PNG)"),
        imageOutput("color_hm_palette", width = "6in", height = "4in")
      ),
      
      tabPanel(
        title = "Clustering",
        value = "Clustering",
        fluidRow( #make each plot smaller to fit more
          p("This tab displays the activity of selected transcription factors"),
          p("- Cells are plotted in 2D according to UMAP dimensionality reduction algorithm"),
          p("- The top row is colored by cluster and the bottom row is colored by transcription factor activity"),
          p("- Only the first two transcriptions factors are displayed"),
          column(width = 10, plotOutput("color_by_cluster", width = "6in", height = "7in"))
        ),
        fluidRow(
  
          column(width = 6, plotOutput("cluster1",width = "4.2in", height = "4in"),
                 downloadButton("download_UMAP_1", "Transcription Factor 1 Activity Plot (PNG)")),

          column(width = 6, plotOutput("cluster2", width = "4.2in",height = "4in"),
                 downloadButton("download_UMAP_2", "Transcription Factor 2 Activity Plot (PNG)")),

        )
      ),
      # tabPanel("Heatmap and Clustering",
      # 
      #          plotOutput("heatmap_cell"),
      #          downloadButton("download_hm_cell", "Heatmap by cell (Png)"),
      #          plotOutput("heatmap_cluster"),
      #          downloadButton("download_hm_cluster", "Heatmap by cluster (Png)"),
      #          imageOutput("color_hm_palette", width = "6in", height = "4in"),
      # 
      #          fluidRow(
      #           textOutput("cluster_UMAP_desc"),
      #           column(width = 8, plotOutput("cluster1",width = "5in", height = "5in"),
      #                  downloadButton("download_UMAP_1", "UMAP scatterplot 1 (Png)")),
      # 
      #           column(width = 8, plotOutput("cluster2", width = "5in",height = "5in"),
      #                  downloadButton("download_UMAP_2", "UMAP scatterplot 2 (Png)")),
      # 
      #          ),
      # 
      #          value = "Heatmap and Clustering"
      # ),
      tabPanel("Time Series",
               textOutput("tf_timeseries_desc"),
               textOutput("timeseries_desc"),
               
               fluidRow(
               plotlyOutput("timeseries1"),
               downloadButton("download_ribbon_1", "Timeseries ribbon plot (Png)"),
               plotOutput("timeseries2"),
               imageOutput("timeseries_color"),
               #plotOutput("timeseries3"),
               #plotlyOutput("timeseries4")
               ),
               value = "Time Series"),
      id = "tabs"
    ))
  ),
  
  
)