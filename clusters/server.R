
# For the plot cache
# NOTE TO SELF:
# To get this to work, need to:
# 1) Make the cache folder, i.e. mkdir cache
# 2) Give the shiny user permissions, i.e. chmod -R a=rwx cache
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
      "gene" = input$gene
    )
    
  })
  
  bubble_input <- reactive({
    
    req(input_dendrogram())
    prep_bubbleplot_input(gene = input_dendrogram()$gene)
    
  })
  
  output$bubble <- renderPlot({
    
    req(bubble_input())
    
    bubbleplot(df = bubble_input())
    
  },
  width = 1171,
  
  # Customize the height of the bubbleplot based on the number of genes which
  # are being displayed, after allocating a baseline height for the x-axis
  height = function() 100 + 30 * length(input_dendrogram()$gene))
  
  output$bubble_hover_info <- renderUI({
    
    # This mouseover tooltip is created following this example
    # https://gitlab.com/snippets/16220
    
    hover <- input$bubble_hover
    
    point <- nearPoints(bubble_input(),
                        hover,
                        xvar = "Cluster",
                        yvar = "Gene_padded",
                        maxpoints = 1) %>% 
      select(Gene, Sample, Cluster, Cell_type, Cell_class, N_cells,
             Expression, Pct1)
    
    # Don't show the tooltip if mouse is not hovering over a point
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
                    "left:", left_px + 0.5, "px; top:", top_px + 0.5, "px; width: 350px;")
    
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
      "gene"   = input$gene_region,
      "region" = input$region
    )
    
  })
  
  # Generate bubble plot and save the output so that we can later split
  # into the plot itself, and the legend
  ribbon <- reactive({
    
    req(input_timecourse())
    
    ribbon_plot(gene   = input_timecourse()$gene,
                region = input_timecourse()$region)
    
  })
  
  # Here, we're grabbing only the plot part, not the legend
  output$plotRibbon <- renderPlot({ ribbon() +
      theme(legend.position = "none")
  })
  
  # Extract the legend to plot separately
  output$ribbonLegend <- renderPlot({
    
    leg <- cowplot::get_legend(ribbon())
    plot_grid(leg)
    
  })
  
}

