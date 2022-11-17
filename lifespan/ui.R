library(shiny)
library(shinycssloaders)

source("../www/ui_functions.R")

# Data directory for local development, ** should be commented on the server **
# data_dir <- "/Volumes/My Passport/braindex/brainspan/"

data_dir <- "."

load(file.path(data_dir, "data/geneName.Rda"))

shinyUI(bootstrapPage(
  
  # Custom styling
  includeCSS("../www/minimal.css"),
  
  navigation(),
  
  beginPage(),
  
  # Application title
  titlePanel("Gene expression in the brain across the lifespan",
             windowTitle = "Lifespan"),
  
  
  # Sidebar with input
  sidebarLayout(
    sidebarPanel(width = 3,
                 selectInput("gene", "Gene", choices = geneName, multiple = FALSE)
    ),
    
    # Output plots
    mainPanel(tabsetPanel(
      
      tabPanel("BrainSpan atlas",

	       tags$br(),

	       p("Visualize the expression of a gene of interest across the lifespan in samples from the BrainSpan project"),
	       
	       p("This analysis is based on code for curve fitting and visualization from Marie Forest and Claudia Kleinman"),

	       p("• In the top row, each curve corresponds to the smoothed expression for one brain region"),

	       p("• Below, one plot is displayed per region, with points representing individual samples -- grey ppints represent observations, and red points indicate imputed values"),
	       p("• The x-axis is log (post-conception weeks), but the the labels correspond to post-conception weeks prenatally, and years postnatally"),

               ws(plotOutput("allcurves", width = "6in", height = "4.5in")),
               
               # Add some whitespace for more breathing room
               br(),
               br(),
               
               ws(plotOutput("regcurves", width = "7in", height = "8in")),
               
               h3("Citation guidelines"),
	       p("These visualizations are based on work from Marie Forest and Claudia Kleinman, please consult for Claudia Kleinman for citation and credit"),
               p("Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. doi:10.1038/nature13185"),
               p("© 2010 Allen Institute for Brain Science. Allen Human Brain Atlas. Available from: https://www.brainspan.org/")
      )
      
    )
    
    )
  ),
  
  # Custom styling
  endPage()
  
)
)
