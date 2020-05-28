
# For the plot cache
# NOTE:
# To get this to work, need to:
# 1) Make the cache folder, i.e. mkdir cache
# 2) On the server, give the shiny user permissions, i.e. chmod -R a=rwx cache
shinyOptions(cache = diskCache("./cache"))

library(cowplot)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)

source("functions.R")
source("../style.R")

# Set a default ggplot2 theme for the app, from cowplot
ggplot2::theme_set(cowplot::theme_cowplot())

server <- function(input, output, session) {
  
  # Dendrogram tab content ----
  
  # Capture all input from this tab as a list in case we want to add
  # more options in the future
  input_dendrogram <- eventReactive(input$update_dendrogram, {
    
    list(
      "scale"  = input$bubble_scale,
      "size"   = input$bubble_size
    )
    
  })
  
  # Generate the input dataframe for the bubbleplot
  bubble_input <- reactive({
    
    req(input_dendrogram())
    bubble_prep(gene  = input$gene,
                scale = input_dendrogram()$scale)
    
  })
  
  # Generate the bubbleplot
  output$bubble <- renderPlot({
    
    req(bubble_input())
    
    bubble_plot(df = bubble_input(),
               max_point_size = input_dendrogram()$size)
    
  },
  
  # Choose width to align horizontally with dendrogram image
  width = 1175,
  
  # Customize the height of the bubbleplot to scale with the number of genes which
  # are being displayed, after allocating a baseline height for the x-axis & legend
  height = function() 150 + 30 * length(input$gene))
  
  # Create a tooltip with cluster / expression information that appears when
  # hovering over a bubble
  # 
  # This adapted from this example https://gitlab.com/snippets/16220
  output$bubble_hover_info <- renderUI({
    
    hover <- input$bubble_hover
    
    # Find the nearest data point to the mouse hover position
    point <- nearPoints(bubble_input(),
                        hover,
                        xvar = "Cluster",
                        yvar = "Gene_padded",
                        maxpoints = 1) %>% 
      select(Gene, Sample, Cluster, Cell_type, Cell_class, N_cells,
             Expression, Pct1)
    
    # Hide the tooltip if mouse is not hovering over a bubble
    if (nrow(point) == 0) return(NULL)
    
    # Create tooltip for mouseover
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # Create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px, "px; top:", top_px, "px; width: 350px;")
    
    # Actual tooltip created as wellPanel, specify info to display
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Gene: </b>",       point$Gene, "<br/>",
                    "<b> Cluster: </b>",    point$Cluster, "<br/>",
                    "<b> Cell type: </b>",  point$Cell_type, "<br/>",
                    "<b> Sample: </b>",     point$Sample, "<br/>",
                    "<b> Expression: </b>", point$Pct1 * point$N_cells, " ",
                    point$Gene, "+ cells out of ", point$N_cells, " cells in cluster <br/>")))
    )
  })
  
  # Timecourse tab content ----
  input_timecourse <- eventReactive(input$update_timecourse, {
    
    list(
      # "gene"   = input$gene_region,
      "region" = input$region
    )
    
  })
  
  # Generate ribbon plot and save the output so that we can later split
  # into the plot itself, and the legend
  ribbon <- reactive({
    
    req(input_timecourse())
    
    ribbon_plot(gene   = input$gene,
                region = input_timecourse()$region)
    
  })
  
  # Grabbing only the plot part, remove the legend
  output$plotRibbon <- renderPlot({ ribbon() +
      theme(legend.position = "none")
  })
  
  # Extract the ribbon plot legend to plot separately
  output$ribbonLegend <- renderPlot({
    
    leg <- cowplot::get_legend(ribbon())
    plot_grid(leg)
    
  })
  
}

