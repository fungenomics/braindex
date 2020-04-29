library(shiny)

source("functions_fda_custom_plotting.R")
source("functions_fda.R")

shinyServer(function(input, output) {

    fda_calc <- reactive({createFDAcurves(prepGene(input$gene))})

    allcurves_plot <- reactive({

        req(fda_calc())
        plotFDAcurves(fda_calc())

    })

    regcurves_plot <- reactive({

        req(fda_calc())
        plotFittedCurvesPerRegion(fda_calc())

    })

    output$allcurves <- renderPlot({allcurves_plot()})
    output$regcurves <- renderPlot({regcurves_plot()})

})
