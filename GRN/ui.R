source("../www/ui_functions.R")

ui <- fluidPage(
  #introjsUI(),
  #useShinyjs(),
  includeCSS("../www/minimal.css"),
   
  navigation(),
   
  beginPage(),
  
  # Application title
  
  titlePanel("Transcription Factor Activity in Single-cell Developmental Atlas",
             windowTitle = "GRN"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # choose which datasets to analyze for the whole app
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
                      
          checkboxInput(inputId = "label", label = "Label Target Gene Nodes", value = FALSE)
         
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
      conditionalPanel(condition = "input.tabs == 'Clustering'",
                       
                       radioButtons("dim_red", "Dimension Reduction Method",
                                          choices = c("UMAP" = "umap",
                                                      "TSNE" = "tsne",
                                                      "PCA" = "pca"),
                                          selected = "umap"),
                       
                       checkboxInput(inputId = "cluster_label", label = "Show Cluster Labels", value = TRUE)
                       ),
      # 3. time series plot
      conditionalPanel(condition = "input.tabs == 'Time Series'"),
      
      # Update everything
      actionButton("update", label = "Update"),
      bookmarkButton(),
    ),
    mainPanel(
      tabsetPanel(
      
        tabPanel(
          strong("This tab displays information corresponding to the selected transcription factors and their predicted target genes.") %>% p(),
          p("• Strength of Association represents the weight of the putative regulatory links between transcription factor and a gene target, 
            as predicted with Genie3, with a higher value indicating a more likely regulatory link."),
          p("• The number of motifs for each gene is identified via the RcisTarget package. The best motif and its sequence logo is displayed."),
          title = "Transcription Factor Target Information",
          #textOutput("general_desc"),
          dataTableOutput("table1"),
          value = "Transcription Factor Target Information"
        ),
        
        tabPanel(
          strong("This tab displays a network visualisation of the inferred regulatory relationship between transcripton factors and target genes.") %>% p(),
          p("• Transcription factors and target genes are represented as nodes with regulatory relationships represented as edges."),
          p("• User input transcription factors are shown in blue. Target genes that are present in the current
            network can be highlighted in orange based on user selection or based on a file containing a gene list."),
          p("• Click on the \"Label Target Gene Nodes\" option to see the label of every gene target and enable a hover option.
            Currently, hover only displays gene name but more information to come soon!"),
          p("• Transcription factors that self regulate are not displayed (i.e no self loops)."),
          p("• File input format: single column csv file with the first row titled 'Gene' and the remaining rows containing a list of genes of interest."),
          title = "Regulatory Network Visualization",
          #textOutput("general_desc"), # this line breaks things/ probably cause you can't have 2 general_desc
          #textOutput("desc"),
          plotlyOutput("network"),
          br(),#so the plotly doesn't overlap with the download button
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          downloadButton("download_network", "Network Visualisation Download (PDF)"),
          #need to work on visualization with the ggNet package
         
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
        strong("This tab displays a heatmap of user selected transcription factor activity per cluster") %>% p(),
        p("• Joint clusters are clusters classified based on the combined data from every developmental 
          time-point per brain region (forebrain or pons); sample cluster are identified based on data from each 
          individual time-point per brain region."),
        title = "Heatmap",
        value = "Heatmap",
        fluidRow(
          plotOutput("heatmap_joint")
        ),
        downloadButton("download_hm_joint", "Heatmap by Joint Cluster (PDF)"),
        div(style = "margin-left: 1.3em; margin-right: 1.3em;",
        fluidRow(
          plotOutput("heatmap_cluster")
        )),
        downloadButton("download_hm_cluster", "Heatmap by Timepoint Cluster (PDF)"),
        imageOutput("color_hm_palette", width = "6in", height = "4in")
      ),
      
      tabPanel(
        title = "Transcription Factor Activity, by Region",
        value = "Clustering",
        fluidRow( #make each plot smaller to fit more
          p("This tab displays the activity of selected transcription factors") %>% strong(),
          p("• Cells are plotted in 2D according to selected dimensionality reduction algorithm"),
          p("• The top row is colored by joint cluster and the bottom row is colored by transcription factor activity"),
          p("• Only the first two transcriptions factors are displayed"),
          column(width = 10, plotOutput("color_by_cluster", width = "6in", height = "7in"))
        ),
        fluidRow(
  
          column(width = 6, plotOutput("cluster1",width = "4.2in", height = "4in"),
                 downloadButton("download_UMAP_1", "Transcription Factor 1 Activity Plot (PDF)")),

          column(width = 6, plotOutput("cluster2", width = "4.2in",height = "4in"),
                 downloadButton("download_UMAP_2", "Transcription Factor 2 Activity Plot (PDF)")),

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
               p("This plot quantifies the proportion of cells (from 0 to 1) at each timepoint where a given
                 transcriptional factor is active, broken down by cell type, to allow for visualizing activity 
                 across the timecourse") %>% strong(),
               
               textOutput("timeseries_desc"),
               br(),
               textOutput("tf_timeseries_desc"),

               fluidRow(
               plotlyOutput("timeseries1"),
               downloadButton("download_ribbon_1", "Timeseries ribbon plot (PDF)"),
               plotOutput("timeseries2"),
               imageOutput("timeseries_color"),
               #plotOutput("timeseries3"),
               #plotlyOutput("timeseries4")
               ),
               value = "Time Series"),
      id = "tabs"
    ))
  ),
  
  # Custom styling
  endPage()
)