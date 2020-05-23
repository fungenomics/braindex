library(shiny)

source("functions_fda_custom_plotting.R")
source("functions_fda.R")
source("../style.R")

shinyServer(function(input, output) {
  
  # First reactive expression based on the input
  # This command runs functional data analysis for a given gene, to find
  # smoothed curves
  fda_calc <- reactive({ createFDAcurves(prepGene(input$gene)) })
  
  # Other reactives depend on first
  # Plot the curves across brain regions in one plot
  allcurves_plot <- reactive({
    
    req(fda_calc())
    
    
    plotFDAcurves(fda_calc())
    
  })
  
  # Plot the curves in each brain region separately, with points
  regcurves_plot <- reactive({
    
    req(fda_calc())
    plotFittedCurvesPerRegion(fda_calc())
    
  })
  
  output$allcurves <- renderPlot({allcurves_plot()})
  output$regcurves <- renderPlot({regcurves_plot()})
  
})
