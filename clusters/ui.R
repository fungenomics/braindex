
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
                 conditionalPanel(condition = "input.tabs == 'joint' || input.tabs == 'sample'",
                                  
                                  selectInput("dr", "Dimensionality reduction",
                                              multiple = FALSE,
                                              choices = c("tSNE", "PCA", "UMAP"),
                                              selected = "tSNE"),
                                  
                                  selectInput("label_clusters", "Label clusters",
                                              choices = c(TRUE, FALSE),
                                              selected = "FALSE"),
                                  
                                  selectInput("feature_palette", "Colour palette",
                                              choices = list("Grey-red"   = "redgrey",
                                                             "Blue-red"   = "rdbu",
                                                             "Blues"      = "blues"),
                                              selected = "redgrey"),
                                  
                                  selectInput("vln_joint_points", "Show points in violin plot",
                                              choices = c(TRUE, FALSE),
                                              selected = FALSE)
                                  
                 ),
                 
                 conditionalPanel(condition = "input.tabs == 'joint'",
                                  
                                  selectInput("dr_clustering", "Annotate cells by",
                                              multiple = FALSE,
                                              choices = c(
                                                "Clustering at the region level" = "joint",
                                                "Clustering at the sample level" = "sample",
                                                "Timepoint"                      = "timepoint"
                                              ),
                                              selected = "Clustering at the region level")
                                  
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
                                       hover = hoverOpts(id = "bubble_hover", clip = FALSE))),
                   
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
               
               fluidRow(
                 
                 plotOutput("scatter_joint", width = "10in", height = "4in")
                 # plotOutput("dr_joint", width = "4.5in", height = "4in",
                 #            hover = hoverOpts(id = "dr_joint_hover", clip = TRUE)),
                 # 
                 # uiOutput("dr_joint_hover_info"),
                 # 
                 # plotOutput("feature_joint", width = "5.33in", height = "4in",
                 #            hover = hoverOpts(id = "feature_joint_hover", clip = TRUE)),
                 # 
                 # uiOutput("feature_joint_hover_info")
                 
                 
               ),
               
               # fluidRow(
               #   
               # 
               #   
               # ),
               
               fluidRow(
                 
                 plotOutput("vln_joint", width = "11in", height = "4in")
                 
               ),
               
               # Specify the value to use when checking if this tab is selected
               value = "joint"
      ),
      
      tabPanel("Single-cell expression, by sample",
               
               p("Single-sample"),
               
               # fluidRow(
               #   
               #   plotOutput("dr_sample", width = "15in", height = "3in")
               #   
               # ),
               # 
               # fluidRow(
               #   
               #   plotOutput("feature_sample", width = "15in", height = "3in")
               #   
               # ),
               # 
               # fluidRow(
               #   
               #   plotOutput("vln_sample", width = "15in", height = "3in")
               #   
               # ),
               
               value = "sample"
               
      ),
      
      id = "tabs"
      
    ))),
  
  # Custom styling
  endPage()
  
)
