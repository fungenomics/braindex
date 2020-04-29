library(shiny)

# data_dir <- "/Volumes/My Passport/braindex/brainspan/"
data_dir <- "."
load(file.path(data_dir, "data/geneName.Rda"))

shinyUI(fluidPage(

    # Application title
    titlePanel("BrainSpan"),

    # Sidebar with input
    sidebarLayout(
        sidebarPanel(width = 3,
                     selectInput("gene", "Gene", choices = geneName, multiple = FALSE)
        ),

        # Output plots
        mainPanel(tabsetPanel(

            tabPanel("Across regions",
                     plotOutput("allcurves", width = "7in", height = "5in")
            ),

            tabPanel("Per region",
                     plotOutput("regcurves", width = "9in", height = "8in")
            )

        )

        )
    )
)
)
