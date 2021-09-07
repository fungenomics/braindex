
source("../www/ui_functions.R")

ui <- function(request){

  bootstrapPage(
  
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
                 
                 conditionalPanel(condition = '!input.upload',
                                  # Gene input field, shared across tabs
                                  selectizeInput(inputId = "gene", label = "Gene", choices = NULL,
                                              multiple = TRUE)),
                 
                 conditionalPanel(condition = 'input.upload',
                                  # Gene list input with a file, shared across tabs
                                  fileInput(inputId = "genelist", label = "Gene list (.txt, .csv, or .tsv)",
                                            buttonLabel = "Browse...",
                                            multiple = FALSE, 
                                            accept = c(".txt", ".csv", ".tsv"),
                                            placeholder = "No file selected")),
  
                 materialSwitch("upload", "Use gene list from file",
                                # status doesn't have any effect other than color scheme. See bootstrap status values
                                status = "success", 
                                value = FALSE,
                                right = TRUE),

                 # Input for dendrogram tab, expression table, and ranked clusters tab
                 conditionalPanel(condition = "(input.tabs == 'dendrogram' || input.tabs == 'exp_table' || input.tabs == 'rank_exp') &&
                                  (input.gene.length > 1 || input.upload)",
                                  materialSwitch("mean_exp", "Display mean expression over the selected genes",
                                                 # status doesn't have any effect other than color scheme. See bootstrap status values
                                                 status = "success",
                                                 value = FALSE,
                                                 right = TRUE),
                 ),
                 
                 # Input for dendrogram tab
                 conditionalPanel(condition = "input.tabs == 'dendrogram'",
                                  
                                  selectInput("bubble_scale", "Scaling",
                                              choices = c("Scale each gene to [0, 1]" = TRUE,
                                                          "Conserve scale across genes" = FALSE),
                                              selected = "Scale each gene to the same range, [0, 1]"),
                                  
                                  sliderInput("bubble_size", "Max point size", min = 3, max = 6,
                                              step = 0.5, value = 4)
                                  
                 ),
                 

                 # Input for all tabs other than dendrogram, expression table, ranked clusters, and heatmap
                 conditionalPanel(condition = "input.tabs != 'dendrogram' && input.tabs != 'exp_table'
                                  && input.tabs != 'rank_exp' && input.tabs != 'heatmap'",
                                  
                                  # Specify the visible label as well as the internal
                                  # strings used to refer to each region, matching
                                  # the paths/files under the data directory
                                  radioButtons("region", "Brain region",
                                               choices = c("Forebrain" = "joint_cortex",
                                                           "Pons" = "joint_pons"))
                 ),
                 
                 # Input for heatmap tab
                 conditionalPanel(condition = "input.tabs == 'heatmap'",
                                  checkboxGroupInput("heatmap_cells", label = "Cell type(s)", 
                                                     choices = list("Progenitors/cycling" = "Progenitors/cyc.",
                                                                    "Oligodendrocytes",
                                                                    "Ependymal", 
                                                                    "Astrocytes",
                                                                    "Neurons",
                                                                    "Non-neuroectoderm" = "Non-neuroect.",
                                                                    "Other"),
                                                     selected = c("Oligodendrocytes",
                                                                  "Ependymal", 
                                                                  "Astrocytes",
                                                                  "Neurons")),
                                  checkboxGroupInput("heatmap_anno", label = "Column annotation(s)", 
                                                     choices = list("Cell class" = "Cell_class",
                                                                    "Brain region" = "Region",
                                                                    "Time point" = "Timepoint"
                                                                    ),
                                                     selected = c("Cell_class"))
                                  ),
                 
                 # Input for tabs on joint analysis by region or by sample
                 conditionalPanel(condition = "input.tabs == 'joint' || input.tabs == 'sample'",
                                  
                                  materialSwitch("vln_points", "Show points in violin plots",
                                                 # status doesn't have any effect other than color scheme. See bootstrap status values
                                                 status = "success", 
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
                                  h4("Dimensionality reduction plots", 
                                            style = "font-size:16px;"),
                                  br(),
                                  
                                  materialSwitch("label_clusters", "Label clusters",
                                                 # status doesn't have any effect other than color scheme. See bootstrap status values
                                                 status = "success", 
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
                 
                 # Update button for all sidebar inputs. Coloured to differentiate
                 # from the bookmark button beside it
                 # tags$head(
                 #   tags$style(HTML('#update{background-color:#4863A0; 
                 #                   color:#FFFFFF;}'))
                 # ),
                 actionButton("update", label = "Update"),
                 
                 # Bookmark button to store current inputs / state of app
                 bookmarkButton()
                 
    ),
    
    # Output plots and tables
    mainPanel(tabsetPanel(
      
      #### ---- Dendrogram tab output ---- 
          
      tabPanel("Dendrogram",
               
               tags$br(),
               tags$b("This tab displays the mean expression of up to 20 genes over each cluster in the mouse scRNAseq development atlas."),
               tags$br(),
               tags$br(),
               p("• Clusters are ordered according to the dendrogram which represents a molecular taxonomy of all cell populations"),
               
               p("• Below the dendrogram, clusters are annotated by brain region, time point, and a cell cycle G2/M phase score"),
               
               p("• Bubble colour encodes the mean expression within the cluster, and bubble size encodes the proportion of cells within each cluster that express the gene"),
               
               p("• Hover over each bubble, or move to the tab containing the table, to get additional details about each cluster & its expression level"),
               
               p("• When selecting more than one gene, use the sidebar switch to plot the mean expression over these genes in a new row of the bubble plot. Note that Pct values are disregarded here, so all bubbles in this row are the same size"),
               
               # Display the image of the cluster dendrogram as in Fig 1 of Jessa et al,
               # Nat Genet, 2019
               div(style = "margin-top: 3em; margin-bottom: -2em !important;",
                   fluidRow(tags$img(src = "tree.png", width = "1150", height = "163"))
               ), 
               
               # Display the bubbleplot
               div(style = "margin-top: 2em; margin-left: 1em; margin-bottom: -5em;
                   overflow-x: visible; overflow-y: visible;",
                   
                   fluidRow(
                     # Set cellWidths equal to the actual width of each plot (server.R)
                     splitLayout(cellWidths = c(1103, 200),
                       
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
               HTML("<br><br><br>"), 
               # Specify the value to use when checking if this tab is selected
               value = "dendrogram"
               
      ),
        
      #### ---- Expression table tab output ---- 
      
      tabPanel("Expression table", 
               
               tags$br(),
               tags$b("This table compares the expression of up to 20 genes in each cluster from the mouse scRNAseq development atlas."),
               tags$br(),
               tags$br(),
               p("• The value in each gene column denotes the mean gene expression per cell in the specified cluster (mean expression)"),
               
               p("• When selecting more than one gene, use the sidebar switch to display the mean expression over these genes in a new column of the table"),
               
               p("• Use the download button below the table to obtain a TSV file with mean expression as well as percent cluster expression values"),
               
               #div(style = "overflow-x: scroll; overflow-y: visible;",
                   fluidRow(
                     reactableOutput("cluster_table", width = 1100) %>% ws
                   ),
               #),
               
               HTML("<br>"),
               
               # Only display download button if update has been pressed at least once
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_bubble", 
                                                 "Download data (TSV)"))
               ),
               HTML("<br><br><br>"), 
               
               # Specify the value to use when checking if this tab is selected
               value = "exp_table"
                   
      ),
      
      #### ---- Timecourse tab output ---- 
      
      tabPanel("Timecourse",
               
               tags$br(),
               
               tags$b("This plot quantifies the proportion of cells (from 0 to 1) at each timepoint where a given gene is detected, broken down by cell type, to allow for visualizing expression across the timecourse."),
               tags$br(),
               tags$br(),
               p("• Use the side bar to select which brain region to interrogate"),
               
               p("• Use the switch above the plot to toggle between static and interactive plots (update button not required)"),
               
               p("• As only one gene can be plotted at a time, use the dropdown tool above the plots to choose which of the input genes to display (update button not required)"),
               
               p("• Download the static version of the plot as a pdf using the button below the plot"),
               
               p("• Be aware of the y-axis, which is computed as the max for each gene"),
               
               tags$br(),
               
               fluidRow(
                 column(6, 
                        wellPanel(
                          selectInput("pick_timecourse", "Select gene to display",
                                      c("Please enter a gene"))
                        )
                 ),
                 column(6,
                        wellPanel(
                          materialSwitch("plotly_ribbon", strong("Interactive ribbon plot"),
                                         # status doesn't have any effect other than color scheme. See bootstrap status values
                                         status = "warning", 
                                         value = FALSE, 
                                         right = TRUE),
                        )
                 )
               ),
               
               # Plot a ribbon plot, showing the proportion of cells in which
               # each gene is detected, broken down by cell type, & across
               # the time course. Either interactive (plotly) or static (ggpplot)
               conditionalPanel(condition = "input.plotly_ribbon",
                                
                                # Plot the ribbon plot & legend as a plotly (interactive) plot
                                plotlyOutput("plotlyRibbon", height = "5in", width = "11.5in") %>% ws
                                ),
               
               conditionalPanel(condition = "!(input.plotly_ribbon)",
                                
                                # Plot the ribbon plot & legend as static plots with ggplot2
                                plotOutput("plotRibbon", height = "8.5in", width = "8in") %>% ws
                                ),
               
               # Only display download button if update has been pressed at least once
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_ribbon",
                                                 "Download ribbon plot (PDF)")
                                )
               ),
               HTML("<br><br><br>"), 
               # Specify the value to use when checking if this tab is selected
               value = "timecourse"
      ),
      
      #### ---- Region joint expression tab output ----
      
      tabPanel("Single-cell expression, by region",
               
               tags$br(),
               
               tags$b("Use this tab to explore the expression of one or more genes at the single-cell level per brain region."),
               tags$br(),
               tags$br(),
               p("• In the top row, the cells are plot in 2D according to a dimensionality reduction algorithm, coloured by cluster (left) or expression (right)"),
               
               p("• If using tSNE or UMAP reduction, hover over the plot coloured by cluster (top left) to identify each cluster. Hover will be disabled if clusters are labeled"),
               
               p("• In the bottom row, violin plots display expression in each cluster, ordered by mean expression"),
               
               p("• If more than one gene is provided, the mean expression of all genes is automatically computed and displayed"),
               
               fluidRow(
                 # plotOutput("scatter_joint", width = "10in", height = "4in") %>% ws
                 
                 splitLayout(cellWidths = c(432, 512), # 432 = 4.5in, 512px = 5.33in
                             #cellArgs = list(style = "padding: 6px"),
                             
                             (plotOutput("dr_joint", 
                                         #width = "4.5in", 
                                         height = "4in",
                                         hover = hoverOpts(id = "dr_joint_hover", clip = TRUE)) %>% ws),
                             
                              (plotOutput("feature_joint", 
                                         #width = "5.33in", 
                                         height = "4in"
                                         #, hover = hoverOpts(id = "feature_joint_hover", clip = TRUE)
                              ) %>% ws)
                         )
                ),
                 
                # Show hover tooltip on clusters if it's not a PCA plot and clusters are unlabeled
                conditionalPanel(condition = "input.dr!='PCA' && !(input.label_clusters)",
                                 fluidRow(uiOutput("dr_joint_hover_info"))
                                ),

                #fluidRow(uiOutput("feature_joint_hover_info")),
               
                fluidRow(
                  plotOutput("vln_joint", width = "11in", height = "4in") %>% ws
                ),
                HTML("<br><br><br>"), 
                # Specify the value to use when checking if this tab is selected
                value = "joint"
      ),
      
      #### ---- Sample joint expression tab output ---- 
      
      tabPanel("Single-cell expression, by sample",
               
               tags$br(),
               
               tags$b("Use these tabs to explore the expression of one or more genes at the single-cell level in each sample."),
               tags$br(),
               tags$br(),
               p("• In the top row of tSNE plots, the cells are plot in the 2D tSNE space and coloured by cluster. In the bottom row, they are instead coloured by expression"),
               
               p("• Each violin plot is coloured by cluster and ordered by the expression level within the given sample"),
               
               p("• If more than one gene is provided, the mean expression of all genes is automatically computed and displayed in both tabs"),
               
               tabsetPanel(
                 
                     tabPanel("tSNE plots",
                         
                         tags$br(),
                         
                         fluidRow(
                           plotOutput("dr_sample", width = "12.5in", height = "2.6in") %>% ws
                         ),
                         
                         fluidRow(
                           plotOutput("feature_sample", width = "12.5in", height = "3in") %>% ws
                         )
                  
                     ),
                
                      tabPanel("Violin plots",
                               
                         tags$br(),
                               
                         fluidRow(
                            plotOutput("vln_sample", width = "10in", height = "20in") %>% ws
                         )
                      )
                 
        ),
        HTML("<br><br><br>"), 
        # Specify the value to use when checking if this tab is selected       
        value = "sample"
      ),
      
      #### ---- Clusters ranked by expression tab output ---- 
      
      tabPanel("Clusters ranked by expression",
               
               tags$br(),
               
               tags$b("This plot displays the mean expression of the selected gene in each cluster, ranked from highest to lowest expression."),
               tags$br(),
               tags$br(),
               p("• The ticks below the plot x-axis provide a general categorization by cell type"),
               
               p("• Use the dropdown tool above the plot to choose which of the input genes to display (update button not required), or use the sidebar toggle to display mean expression over all input genes (update button required)"),
               
               p("• Be aware of the y-axis, which is bounded by the maximum expression value present"),
               
               p("• Use the horizontal scroll bar at the bottom of the plot to view the plot's full width"),
               
               p("• Download the plot as a PDF file using the button below the plot"),
               
               HTML("<br>"),
               
               fluidRow(
                 column(6, 
                        conditionalPanel(condition="!input.mean_exp",
                                        wellPanel(
                                          selectInput("pick_ranked", "Select gene to display",
                                                      c("Please enter a gene"))
                                        )
                        )
                 )
               ),
               
               # Plot will allow scrolling to view the full horizontal width of the plot
               div(style = "width: 1152px; overflow-x: auto; overflow-y: visible;",
                 fluidRow(
                   plotOutput("rank_tick_plot", width = "17in", height = "5in") %>% ws
                 ) 
               ),
               
               HTML("<br>"),
               
               # Only display download button if update has been pressed at least once
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_ranked_plot",
                                                 "Download ranked plot (PDF)")
                                )
               ),
               
               HTML("<br><br><br>"), 
               # Specify the value to use when checking if this tab is selected       
               value = "rank_exp"
      ),
      
      #### ---- Cell types clustered by expression tab output ---- 
      
      tabPanel("Cell types clustered by expression",
               
               tags$br(),
               tags$b("This plot is a heatmap clustering input genes and cell types together based on their mean expression within clusters."),
               tags$br(),
               tags$br(),
               p("• The heatmap's hierarchical clustering method requires at least two genes as input. An error message will display if only one gene is provided"),
               
               p("• At least one cell type must be selected: an error message will display if none are checked off"),
               
               p("• The coloured bar above the heatmap provides a categorization of cell clusters by general cell type"),
               
               p("• The tree to the left of the heatmap indicates the clustering of genes, and the tree above the heatmap indicates the clustering of cell types"),
               
               p("• Download the plot as a PDF file using the button below the plot"),
               
               # TODO: implement scroll accommodating variable widths:
               # overflow-x: auto doesn't hide the scroll bar for smaller plots...
               #div(style = "width: 1500px; overflow-x: visible; overflow-y: visible;",
                 fluidRow(
                   uiOutput("heatmapUI") %>% ws
                 ),
               #),
               
               HTML("<br>"),
               
               # Only display download button if update has been pressed at least once
               conditionalPanel(condition='input.update!=0',
                                fluidRow(
                                  downloadButton("download_heatmap",
                                                 "Download heatmap (PDF)")
                                )
               ),
               
               HTML("<br><br><br>"), 
               # Specify the value to use when checking if this tab is selected       
               value = "heatmap"
      ),
      
      id = "tabs"
      
    ))),
  
  # Custom styling
  endPage()
  
)}
