
library(readr)
source("../www/ui_functions.R")

# Load names of genes detected in mouse to provide choices in input
genes_mouse <- data.table::fread("data/joint_mouse/joint_mouse.gene_names.tsv", data.table = FALSE)$genes

ui <- bootstrapPage(
  
  # Custom styling
  includeCSS("../www/minimal.css"),
  
  navigation(),
  
  beginPage(),
  
  # Application title
  titlePanel("Expression in single-cell developmental atlases, by cluster",
             windowTitle = "Clusters"),
  
  # Sidebar with input
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 # Gene input field, shared across tabs
                 selectInput("gene", "Gene", choices = genes_mouse,
                             multiple = TRUE),
                 
                 # Input for dendrogram tab
                 conditionalPanel(condition = "input.tabs == 'dendrogram'",
                                  
                                  selectInput("bubble_scale", "Scaling",
                                              choices = c("Scale each gene to [0, 1]" = TRUE,
                                                          "Conserve scale across genes" = FALSE),
                                              selected = "Scale each gene to the same range, [0, 1]"),
                                  
                                  sliderInput("bubble_size", "Max point size", min = 3, max = 6,
                                              step = 0.5, value = 4)
                                  
                 ),
                 
                 # Input for timecourse tab
                 conditionalPanel(condition = "input.tabs != 'dendrogram'",
                                  
                                  # Specify the visible label as well as the internal
                                  # strings used to refer to each region, matching
                                  # the paths/files under the data directory
                                  radioButtons("region", "Brain region",
                                               choices = c("Forebrain" = "joint_cortex",
                                                           "Pons" = "joint_pons"))
                 ),
                 
                 # Input for joint analysis tab
                 conditionalPanel(condition = "input.tabs == 'joint'",
                                  
                                  selectInput("dr", "Dimensionality reduction",
                                              multiple = FALSE,
                                              choices = c("tSNE", "PCA", "UMAP"),
                                              selected = "tSNE"),
                                  
                                  selectInput("dr_clustering", "Group cells by",
                                              multiple = FALSE,
                                              choices = c(
                                                "Clustering at the region level",
                                                "Clustering at the sample level"
                                              ),
                                              selected = "Clustering at the region level"),
                                  
                                  selectInput("label_clusters", "Label clusters",
                                             choices = c(TRUE, FALSE),
                                             selected = "FALSE")
                                  
                 ),
                 
                 actionButton("update", label = "Update")
                 
    ),
    
    # Output plots
    mainPanel(tabsetPanel(
      
      tabPanel("Dendrogram",
               
               tags$br(),
               
               p("Display up to 6 genes"),
               
               # Display the image of the cluster dendrogram as in Fig 1 of Jessa et al,
               # Nat Genet, 2019
               div(style = "margin-top: 3em; margin-bottom: -2em !important;",
                   fluidRow(tags$img(src = "tree.png", width = "1150", height = "163"))
               ), 
               
               # Display the bubbleplot
               div(style = "margin-top: 2em; margin-left: 1.3em; margin-bottom: -5em;",
                   fluidRow(plotOutput("bubble",
                                       hover = hoverOpts(id = "bubble_hover", clip = TRUE))),
                   
                   # UI for tooltip
                   fluidRow(
                     uiOutput("bubble_hover_info"))
                   
               ),
               
               fluidRow(DT::dataTableOutput("cluster_table", width = 1100)),
               
               # fluidRow(
               #   downloadButton("download_bubble", "Download data (TSV)")),
               
               # Specify the value to use when checking if this tab is selected
               value = "dendrogram"
               
      ),
      
      tabPanel("Timecourse",
               
               # Plot a ribbon plot, showing the proportion of cells in which
               # each gene is detected, broken down by cell type, across
               # the time course
               plotOutput("plotRibbon"),
               plotOutput("ribbonLegend"),
               
               # Specify the value to use when checking if this tab is selected
               value = "timecourse"
      ),
      
      tabPanel("Single-cell expression, by region",
               
               tags$br(),
               
               # Plot a ribbon plot, showing the proportion of cells in which
               # each gene is detected, broken down by cell type, across
               # the time course
               fluidRow(
                 
                 column(4,
                        plotOutput("dr_joint", width = "5.5in", height = "5in")),
                 
                 column(4,
                        p("Test"))
                 
               ),
               
               # Specify the value to use when checking if this tab is selected
               value = "joint"
      ),
      
      id = "tabs"
      
    ))),
  
  # Custom styling
  endPage()
  
)
