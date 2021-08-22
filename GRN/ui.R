source("../www/ui_functions.R")

ui <- fluidPage(
  #introjsUI(),
  useShinyjs(),
  includeCSS("../www/minimal.css"),
   
  navigation(),
   
  beginPage(),
  
  # Application title
  
  titlePanel("Transcription Factor Activity in Single-cell Developmental Atlas",
             windowTitle = "GRN"),
  # ---------------- Side Panel --------------------------------------------- 
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
      selectizeInput(inputId = "TF",
                     label = "Transcription Factor Input",
                     choices = NULL,
                     multiple = TRUE,
                     selected = c("Arx", "Lef1")),
      #text output for user TF input
      htmlOutput("tf_check"),
      br(),

#----------------------------ggNet visualisation---------------------------------------------
      conditionalPanel(
          condition = "input.tabs == 'Regulatory Network Visualization'",
          selectizeInput(inputId = "gene",
                         label = "Genes of Interest",
                         choices = NULL,
                         multiple = TRUE,
                         selected = c("Dlx6","Sox6")),
          #allows user to input a file containing a list of genes to query
          #this input returns a dataframe with 4 columns: name, size, type, datapath
          #each file entered is one row
          
          uiOutput('file1_ui'), ## instead of fileInput('file1', label = NULL) so that the file can be reset 

          actionButton("reset", label = "Reset File"),
                      
          checkboxInput(inputId = "label", label = "Label Target Gene Nodes", value = FALSE)
         
      ),
# -----------------TF info table ---------------------------------------------
      conditionalPanel(
          condition = "input.tabs == 'Transcription Factor Target Information'",
  
      ),
     
# ----------------Heatmap ---------------------------------------------
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
                                   label = "Time-point to Visualize",
                                   choices = c("All Time-Points", dev_time_points),
                                   multiple = FALSE,
                                   selected = "All Time-Points")
                                          
                       ),
# -----------------DR plots ---------------------------------------------
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
#----------------------Active specific-----------------      
      conditionalPanel(condition = "input.tabs == 'active_specific'",
                       #checkboxInput(inputId = "use_in",
                                     #label = "View Selected Regulons"),
                       selectizeInput(inputId = "active_specific_cluster",
                                      label = "Cluster of Interest",
                                      choices = NULL,
                                      multiple = FALSE),
                       sliderInput("fc", "Fold Change Cut-off",
                                   min = 1, max = 4, value = 1.5, step = 0.25, ticks = TRUE),
                       #materialSwitch("dendro", "See Dendrogram", status = "success", right = TRUE)
                       ),
      
      # Update everything
      actionButton("update", label = "Update"),
      
      #save url to server button
      bookmarkButton(),

      
    ),
# -----------------Main Panel ---------------------------------------------
    mainPanel(
      tabsetPanel(
        
        #-----------------------Bubble-plot-------------------
        tabPanel(
         title = "Dendrogram",
         tags$br(),
         tags$b("This tab displays the activity and activity fold change of up to 20 TFs over each cluster in the selected brain region."),
         tags$br(),
         tags$br(),
         p("• Clusters are ordered according to the dendrogram which represents a molecular taxonomy of all cell populations"),
         
         p("• Below the dendrogram, clusters are annotated by brain region, and time point."),
         
         p("• Bubble colour encodes the mean TF activity within the cluster, and bubble size encodes the fold change of TF activity in each cluster compared to all other clusters - effectively describing TF specificity."),
         
         p("• Hover over each bubble to get additional details about each cluster & its expression level"),
         
        uiOutput("dend_image"),
         
         div(style = "margin-top: 2em; margin-left: 0em; margin-bottom: -5em;
                   overflow-x: visible; overflow-y: visible;",
             
             fluidRow(
               # Set cellWidths equal to the actual width of each plot (server.R)
               splitLayout(cellWidths = c(1148, 200),
                           
                           # Bubble plot(s)
                           (plotOutput("bubble",
                                       hover = hoverOpts(id = "bubble_hover", clip = FALSE),
                                       height = 2000) %>% ws),
                           
                           # Gene labels
                           # No spinner to prevent confusing user, because there is only 1 plot
                           (plotOutput("bubble_labels", height = 2000)) 
               )
               
             ),
             
             # UI for tooltip
             fluidRow(
               uiOutput("bubble_hover_info")),
             
         ),
         value = "bubble" 
        ),
        # -----------------Heatmap ---------------------------------------------       
        tabPanel(
          strong("This tab displays a heatmap of user selected TF activity per cluster") %>% p(),
          p("• Values in the heatmap represent the mean TF activity per cluster."),
          p("• Joint clusters are classified based on the combined data from every developmental 
          time-point in a brain region (forebrain or pons); sample cluster are identified based on data from each 
          individual time-point per brain region."),
          p("• Use the \"Time-point to Visalise\" option to select which (if not all) time-points
          to visualise in the sample cluster heatmap."),
          
          fluidRow(
            column(width = 3,  materialSwitch(inputId = "heatmap_toggle", 
                                              label = "Explore per Time-Point Heatmap", 
                                              value = FALSE, status = "success")),
            column(width = 3, actionButton("info2", "What Is This?"))
          ),
          # materialSwitch(inputId = "heatmap_toggle", 
          #                label = "Explore per Time-Point Heatmap", 
          #                value = FALSE, status = "success"),
          htmlOutput("hm_data"),
          title = "TF Activity Heatmap",
          value = "Heatmap",
          fluidRow(
            plotOutput("heatmap_joint")
          ),
          
          downloadButton("download_hm_joint", "Heatmap Download (PDF)"),
          
          (div(style='width:800px;overflow-x: scroll;',
               uiOutput("heatmap_cluster"))),
          
          # div(style = "margin-left: 1.3em; margin-right: 1.3em;",
          # fluidRow(
          #   plotOutput("heatmap_cluster")
          # )),
          downloadButton("download_hm_cluster", "Heatmap Download (PDF)"),
          #imageOutput("color_hm_palette", width = "6in", height = "4in")
        ),
        
        # -----------------DR plots ---------------------------------------------     
        tabPanel(
          title = "TF Activity, by Region",
          value = "Clustering",
          
          p("This tab displays the activity of selected transcription factors") %>% strong(),
          p("• Cells are plotted in 2D according to selected dimensionality reduction algorithm"),
          p("• The top row is colored by joint cluster and the bottom row is colored by transcription factor activity"),
          p("• Only the first two transcriptions factors are displayed"),
          fluidRow(
            column(width = 3, materialSwitch(inputId = "cluster_toggle", 
                                             label = "Explore per Time-Point TF Activity",
                                             value = FALSE, status = "success")),
            column(width = 6, actionButton("info3", "What Is This?"))
          ),
          htmlOutput("dr_data"),
          fluidRow(
            column(width = 10, plotOutput("color_by_cluster", width = "6in", height = "7in"))
          ),
          fluidRow(
            
            column(width = 6, plotOutput("cluster1",width = "4.2in", height = "4in"),
                   downloadButton("download_UMAP_1", "Transcription Factor 1 Activity Plot (PDF)")),
            
            column(width = 6, plotOutput("cluster2", width = "4.2in",height = "4in"),
                   downloadButton("download_UMAP_2", "Transcription Factor 2 Activity Plot (PDF)")),
            
          )
        ),
        
        # -----------------ggNet visualisation ---------------------------------------------       
        tabPanel(
          strong("This tab displays a network visualisation of the inferred regulatory relationship between TFs and target genes.") %>% p(),
          p("• TFs and target genes are represented as nodes with regulatory relationships represented as edges."),
          p("• User input TFs are shown in blue. Target genes that are present in the current
            network can be highlighted in orange based on a user-input gene list, either through the
            \"Genes of Interest\"input or through a file input."),
          p("• Click on the \"Label Target Gene Nodes\" option to see the label of every gene target and enable a hover option.
            Currently, hover only displays gene name but more information to come soon!"),
          p("• TFs that self regulate are not displayed (i.e no self loops)."),
          p("• File input format: single column csv file with the first row titled 'Gene' and the remaining rows containing a list of genes of interest."),
          title = "GRN Visualization", 
          
          fluidRow(
            column(width = 3, materialSwitch(inputId = "grn_toggle", label = "Explore per Time-Point GRN",
                         value = FALSE, status = "success")),
            column(width = 3, actionButton("info", "What Is This?"))
            ),
          
          htmlOutput("grn_data"),
          
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
         
          value = "Regulatory Network Visualization"
        ),
        # -----------------TF info table ---------------------------------------------      
        tabPanel(
          strong("This tab displays information corresponding to the selected TFs and their inferred target genes.") %>% p(),
          p("• Strength of Association represents the weight of the putative regulatory links between transcription factor and a gene target, 
            as predicted with Genie3, with a higher value indicating a more likely regulatory link."),
          p("• The number of motifs for each gene is identified via the RcisTarget package. The best motif and its sequence logo is displayed."),
          title = "TF Target Information",
          
          fluidRow(
            column(width = 3,  materialSwitch(inputId = "table_toggle", label = "Explore per Time-Point Data",
                                              value = FALSE, status = "success")),
            column(width = 3, actionButton("info1", "What Is This?"))
          ),
          # materialSwitch(inputId = "table_toggle", label = "Explore per Time-Point Data",
          #                value = FALSE, status = "success"),
          htmlOutput("table_data"),
          #textOutput("general_desc"),
          dataTableOutput("table1"),
          value = "Transcription Factor Target Information"
        ),

      
   
      # -----------------Time Series ---------------------------------------------
      tabPanel("Time Series",
               p("This plot quantifies the proportion of cells (from 0 to 1) at each timepoint where a given
                 TF is active, broken down by cell type, to allow for visualizing activity 
                 across time.") %>% strong(),
               p("• For any given cell, any given TF is considered active if its activity in that cell
                  is higher than a TF activity threshold determined in the SCENIC pipeline."),
               p("• The time series for the first TF selected in the sidebar will be an interactive plot, with
                 the remaining plots being static."),
               
               textOutput("timeseries_desc"),
               br(),
               #textOutput("tf_timeseries_desc"),

               fluidRow(
               plotlyOutput("timeseries1"),
               downloadButton("download_ribbon_1", "Timeseries ribbon plot (PDF)"),
               plotOutput("timeseries2"),
               #imageOutput("timeseries_color"),
               #plotOutput("timeseries3"),
               #plotlyOutput("timeseries4")
               plotlyOutput("cell_proportion_timeseries")
               ),
               #imageOutput("proportion_timeseries", width = "auto", height = "auto"),
               value = "Time Series"),
      #-----------------------------Active specific------------------------
      tabPanel("TF Activity and Specificity",
               p("This tab displays information on where a TF is the most active and 
                 specific.") %>% strong(),
               p("• The \"By Cluster\" sub-tab visualises the AUC of each TF in the selected cluster on the y-axis
                  and the average AUC in all other clusters in the sample on the x-axis. User-selected
                 fold change cutoff is displayed as a diagonal line and TFs with fold change greater than
                 this cutoff is summarized in the table." ),
               p("• The \"By TF\" sub-tab visualises the AUC of the first 4 selected TFs across all clusters 
                 in the sample. If there are more than 30 clusters in the sample, the 30 clusters with the highest
                 AUC value for that TF is shown."),
               
               fluidRow(
                 column(width = 3, materialSwitch(inputId = "as_toggle", label = "Explore per Time-Point Data",
                                                  value = FALSE, status = "success")),
                 column(width = 3, actionButton("info4", "What Is This?"))
               ),
               htmlOutput("as_data"),
               tabsetPanel(
                 tabPanel("By Cluster", uiOutput("as_clust"), 
                          downloadButton("as_scatter_download", "Scatter-Plot Download (PDF)"), 
                          br(), 
                          fluidRow(
                            column(12, align="center",
                                   plotOutput("active_specific_dr",  width = "4in", height = "4in"),
                            )
                          ),
                          value = "by_clust"),
                 tabPanel("By TF", uiOutput("as_tf"), value = "by_tf"),
                 id = "as_tabs"
               ),
               #uiOutput("as_plots"),
               # fluidRow(
               #   column(width = 7, plotOutput("active_specific_scatter")),
               #   column(width = 5, tableOutput("active_specific_table"))
               #   
               # ),
               plotOutput("active_specific_bar"),
               downloadButton("as_bar_download", "Bar-Plot Download (PDF)"),
               value = "active_specific"),
      id = "tabs"
    ))
  ),
  
  # Custom styling
  endPage()
)