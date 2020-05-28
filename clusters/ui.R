
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
                 
                 # Input for dendrogram tab
                 conditionalPanel(condition = "input.tabs == 'dendrogram'",
                                  
                                  selectInput("gene", "Gene", choices = genes_mouse,
                                              multiple = TRUE),
                                  
                                  actionButton("update_dendrogram", label = "Update")
                 ),
                 
                 # Input for timecourse tab
                 conditionalPanel(condition = "input.tabs == 'timecourse'",
                                  
                                  # TODO: Dynamically provide the right gene
                                  # names as choices based on which brain region
                                  # is selected
                                  selectInput("gene_region", "Gene",
                                              choices = genes_mouse,
                                              multiple = FALSE),
                                  
                                  # Specify the visible label as well as the internal
                                  # strings used to refer to each region, matching
                                  # the paths/files under the data directory
                                  radioButtons("region", "Brain region",
                                               choices = c("Forebrain" = "joint_cortex",
                                                           "Pons" = "joint_pons"),
                                               selected = "Forebrain"),
                                  
                                  actionButton("update_timecourse", label = "Update")
                                  
                 )
                 
    ),
    
    # Output plots
    mainPanel(tabsetPanel(
      
      tabPanel("Dendrogram",
               
               # Display the image of the cluster dendrogram as in Fig 1 of Jessa et al,
               # Nat Genet, 2019
               div(style = "margin-top: 3em; margin-bottom: -2em !important;",
                   fluidRow(tags$img(src = "tree.png", width = "1150", height = "163"))
               ), 
               
               div(style = "margin-top: 3em; margin-left: 1.3em;",
                   fluidRow(plotOutput("bubble",
                                       hover = hoverOpts(id = "bubble_hover", clip = TRUE))),
                   fluidRow(
                     uiOutput("bubble_hover_info"))
                   
               ),
               
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
      
      id = "tabs"
      
    ))),
  
  # Custom styling
  endPage()
  
)
