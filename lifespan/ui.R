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
               ws(plotOutput("allcurves", width = "6in", height = "4.5in")),
               
               # Add some whitespace for more breathing room
               br(),
               br(),
               
               ws(plotOutput("regcurves", width = "7in", height = "8in")),
               
               h3("Citation guidelines"),
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
