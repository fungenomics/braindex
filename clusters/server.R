#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# shinyOptions(cache = diskCache("./cache"))

library(cowplot)
library(shiny)
library(glue)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(grid)
library(gridExtra)

source("shiny_functions.R")
theme_set(theme_cowplot())

load("data/shiny_precomputed_resources.Rda")

genes <- readr::read_tsv("data/joint_cortex.gene_names.tsv")
genes <- genes$genes

server <- function(input, output, session) {
  
  # # Note that the "gene" is really an arbitrary feature including any metadata column
  # updateSelectInput(session = session, "gene", choices =
  #                     list(
  #                       "Quality control stats" = qc_stats,
  #                       "Gene signatures" = gene_sigs,
  #                       "Genes" = all_genes
  #                     )
  # )
  # 
  # observe({
  #   
  #   req(input$gene)
  #   gene <- input$gene
  #   
  #   if (gene %in% gene_sigs) {
  #     
  #     updateSelectInput(session, "gene_stat", selected = "ssgsea")  
  #     updateSelectInput(session, "palette", selected = "rdbu")
  #     
  #   } else if (gene %in% qc_stats) {
  #     
  #     updateSelectInput(session, "gene_stat", selected = "qcstat")  
  #     updateSelectInput(session, "palette", selected = "viridis")
  #     
  #   } else {
  #     
  #     updateSelectInput(session, "gene_stat", selected = "indiv")  
  #     updateSelectInput(session, "palette", selected = "redgrey")
  #     
  #   }
  #   
  # })
  
  # Reactive expressions ----
  
  
  # Expression
  time_plot <- reactive({
    
    req(input$gene, input$mode, input$points, input$min_prop, input$span, input$trend, input$x_axis, input$jitter_width)
    
    if (input$mode == "age") {

    plot_gene_in_time(gene         = input$gene,
                      span         = input$span,
                      points       = input$points,
                      min_prop     = input$min_prop,
                      trend        = input$trend,
                      x_axis       = input$x_axis,
                      jitter_width = input$jitter_width)
   
    } else if (input$mode == "pseudotime") {

	  plot_gene_in_pseudotime(gene         = input$gene,
                      span         = input$span,
                      points       = input$points,
                      min_prop     = input$min_prop,
                      trend        = input$trend)

    } else if (input$mode == "detected") {

	ribbon_plot(gene = input$gene)

    }

  })
  
  # output$time <- renderCachedPlot({
  #   time_plot()$plot
  # },
  # cacheKeyExpr = list(input$gene, input$span, input$min_prop, input$points, input$trend, input$x_axis, input$jitter_width)
  # )
  
  output$time <- renderPlot({time_plot()$plot})

  # output$time_legend <- renderCachedPlot({
   #  leg <- time_plot()$legend
  #   
  #   grid.newpage()
  #   grid.draw(leg)
   #  
  # },
  # cacheKeyExpr = list(input$gene, input$span, input$min_prop, input$points, input$trend, input$x_axis, input$jitter_width)
 #  )
  
	output$time_legend <- renderPlot({

		leg <- time_plot()$legend
		grid.newpage()
		grid.draw(leg)

	})

  # Feature plot download
  output$download_time <- downloadHandler(filename = "time.png",
                                          contentType = "image/png",
                                          content = function(file) {
                                            ggsave(filename = file, plot = time_plot()$plot, width = 10, height = 4)
                                          })
  
}

