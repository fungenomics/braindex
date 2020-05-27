
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
  
  # Display the image of the cluster dendrogram as in Fig 1 of Jessa et al,
  # Nat Genet, 2019
  # output$dendrogram <- renderImage({
  #   list(src = "img/tree.png",
  #        contentType = 'image/png',
  #        width = 965,
  #        height = 200)
  # }, deleteFile = FALSE)
  
  output$bubble <- renderCachedPlot({
    
    req(input_dendrogram())
    
    plot_grid(
      NULL,
      # Generate a bubble plot for expression across clusters in dendrogram order
      bubbleplot_expr(gene = input_dendrogram()$gene),
      # Add a dummy element on the right to customize alignment w/ dendrogram image
      rel_widths = c(0.04, 0.96))
    
  }, cacheKeyExpr = { input_dendrogram()$gene })
  
  # Customize the height of the bubbleplot based on the number of genes which
  # are being displayed, after allocating a baseline height for the x-axis
  # labels
  plotHeight <- reactive(100 + (40 * length(input_dendrogram()$gene)))
  
  # Output element which displays the bubble plot with the reactive height
  output$plotBubble <- renderUI({
    plotOutput("bubble", height = plotHeight(), width = 1199)
  })
  
  
  # Timecourse tab content ---
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

