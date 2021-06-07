
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
  
  #### ---- Sidebar (input) ----
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
                 

                 # Input for all tabs other than dendrogram & table
                 conditionalPanel(condition = "input.tabs != 'dendrogram' && input.tabs != 'exp_table'",
                                  
                                  # Specify the visible label as well as the internal
                                  # strings used to refer to each region, matching
                                  # the paths/files under the data directory
                                  radioButtons("region", "Brain region",
                                               choices = c("Forebrain" = "joint_cortex",
                                                           "Pons" = "joint_pons"))
                 ),
                 
                 # # Input for timecourse ribbon plot tab
                 # conditionalPanel(condition = "input.tabs == 'timecourse'",
                 #                  
                 #                  materialSwitch("plotly_ribbon", "Interactive ribbon plot",
                 #                                 status = "success", # Success status doesn't have any effect other than green color scheme
                 #                                 value = FALSE, 
                 #                                 right = TRUE)
                 #                  ),
                 
                 # Input for tabs on joint analysis by region or by sample
                 conditionalPanel(condition = "input.tabs == 'joint' || input.tabs == 'sample'",
                                  
                                  materialSwitch("vln_points", "Show points in violin plots",
                                                 status = "success", # Success status doesn't have any effect other than green color scheme
                                                 value = FALSE,
                                                 right = TRUE),
                                  
                                  # Input only in region joint analysis tab
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
                                                   
                                  
                                  # Produce separation in sidebar: all options below are about the plots
                                  hr(style = "border-top: 1px solid #000000;"),
                                  h5(strong("Dimensionality reduction plots")),
                                  br(),
                                  
                                  materialSwitch("label_clusters", "Label clusters",
                                                 status = "success",  # Success status doesn't have any effect other than green color scheme
                                                 value = FALSE,
                                                 right = TRUE),
                                  
                                  selectInput("feature_palette", "Expression colour palette",
                                              choices = list("Grey-red"   = "redgrey",
                                                             "Blue-red"   = "rdbu",
                                                             "Blues"      = "blues"),
                                              selected = "redgrey"),
                                  
                                  # Input only in region joint analysis tab
                                  conditionalPanel(condition = "input.tabs == 'joint'",
                                                   
                                                   selectInput("dr", "Dimensionality reduction method",
                                                               multiple = FALSE,
                                                               choices = c("tSNE", "PCA", "UMAP"),
                                                               selected = "tSNE")
                                                   
                                  ),
                 ),
                 
                 actionButton("update", label = "Update")
                 
    ),
    
    # Output plots and tables
    mainPanel(tabsetPanel(
      
      #### ---- Dendrogram tab output ---- 
          
      tabPanel("Dendrogram",
               
               tags$br(),
               p("This tab displays the mean expression of up to 6 genes in each cluster from the mouse scRNAseq development atlas"),
               
               p("• Clusters are ordered according to the dendrogram which represents a molecular taxonomy of all cell populations"),
               
               p("• Below the dendrogram, clusters are annotated by brain region, time point, and a cell cycle G2/M phase score"),
               
               p("• Bubble colour encodes the mean expression, and bubble size encodes the proportion of cells within each cluster"),
               
               p("• Hover over each bubble, or move to the tab containing the table, to get additional details about each cluster & its expression level"),
               
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
               
               # Specify the value to use when checking if this tab is selected
               value = "dendrogram"
               
      ),
        
      #### ---- Expression table tab output ---- 
      
      tabPanel("Expression table", #TODO: confirm a better name
               
               tags$br(),
               p("This table compares the expression of up to 6 genes in each cluster from the mouse scRNAseq development atlas"),
               
               p("• The value in each gene column denotes the mean gene expression per cell in the specified cluster (mean expression)"),
               
               p("• Use the download button below the table to obtain a TSV file with mean expression as well as percent cluster expression values"),
               
               fluidRow(DT::dataTableOutput("cluster_table", width = 1100)),
               
               # Only allow download button to display if update button has been pressed 
               # TODO: figure out how to do this condition in server.R using req() in reactive()
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_bubble", 
                                                 "Download data (TSV)"))
               ),
               
               # Specify the value to use when checking if this tab is selected
               value = "exp_table"
                   
      ),
      
      #### ---- Timecourse tab output ---- 
      
      tabPanel("Timecourse",
               
               tags$br(),
               
               p("This plot quantifies the proportion of cells (from 0 to 1) at each timepoint where a given gene is detected, broken down by cell type, to allow for visualizing expression across the timecourse"),
               
               p("• Use the side bar to select which brain region to interrogate"),
               
               p("• Use the switch below to toggle the plot's interactivity (immediate response)"),
               
               p("• The download button will store the non-interactive version of the plot as a pdf"),
               
               p("• Be aware of the y-axis, which is computed as the max for each gene"),
               
               p("• If more than one gene is provided, only the first gene is plotted"),
               
               materialSwitch("plotly_ribbon", "Interactive ribbon plot",
                              status = "success", # Success status doesn't have any effect other than green color scheme
                              value = FALSE, 
                              right = TRUE
               ),
               
               # Only allow download button to display if update button has been pressed 
               # TODO: figure out how to do this condition in server.R using req() in reactive()
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_ribbon",
                                                 "Download ribbon plot (PDF)")
                                )
               ),
               
               # Plot a ribbon plot, showing the proportion of cells in which
               # each gene is detected, broken down by cell type, across
               # the time course, either interactively or statically
               
               # TODO: make the plot type change only after update button is pressed
               conditionalPanel(condition = "input.plotly_ribbon",
                                
                                # Plot the ribbon plot & legend as a plotly (interactive) plot
                                plotlyOutput("plotlyRibbon", height = "5in", width = "12.5in")
                                
                                ),
               
               conditionalPanel(condition = "!(input.plotly_ribbon)",
                                
                                # Plot the ribbon plot & legend as static plots with ggplot2
                                plotOutput("plotRibbon", height = "10in", width = "10in")
                                
                                ),
               

               
               # Specify the value to use when checking if this tab is selected
               value = "timecourse"
      ),
      
      #### ---- Region joint expression tab output ----
      
      tabPanel("Single-cell expression, by region",
               
               tags$br(),
               
               p("Use this tab to explore the expression of one or more genes at the single-cell level per brain region"),
               
               p("• In the top row, the cells are plot in 2D according to a dimensionality reduction algorithm, coloured by cluster (left) or expression (right)"),
               
               p("• In the bottom row, violin plots display expression in each cluster, ordered by mean expression"),
               
               p("• If more than one gene is provided, the mean expression of all genes is automatically computed and displayed"),
               
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
      
      #### ---- Sample joint expression tab output ---- 
      
      tabPanel("Single-cell expression, by sample",
               
        tabsetPanel(
                 
          tabPanel("tSNE plots",
                   
                   tags$br(),
                   
                   p("Use this tab to explore the expression of one or more genes at the single-cell level in each sample"),
                   
                   p("• In the top row, the cells are plot in the 2D tSNE space, coloured by cluster"),
                   
                   p("• In the bottom row, the cells are plot in the 2D tSNE space, coloured by expression"),
                   
                   p("• If more than one gene is provided, the mean expression of all genes is automatically computed and displayed"),
                   
                   fluidRow(
                     
                     plotOutput("dr_sample", width = "12.5in", height = "2.6in")
                     
                   ),
                   
                   fluidRow(
                     
                     plotOutput("feature_sample", width = "12.5in", height = "3in")
                     
                   )
            
          ),
                
          tabPanel("Violin plots",
                   
                   tags$br(),
                   
                   p("Use this tab to explore the expression of one or more genes at the single-cell level in each sample"),
                   
                   p("• Each violin plot is coloured by cluster and ordered by the expression level within the given sample"),
                   
                   p("• If more than one gene is provided, the mean expression of all genes is automatically computed and displayed"),
                   
                   fluidRow(
                     
                     plotOutput("vln_sample", width = "10in", height = "20in")
                     
                   )
                   
          )
                 
        ),
        
        # Specify the value to use when checking if this tab is selected       
        value = "sample"
      ),
      
      id = "tabs"
      
    ))),
  
  # Custom styling
  endPage()
  
)
