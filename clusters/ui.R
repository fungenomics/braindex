
library(readr)

genes <- readr::read_tsv("data/joint_cortex.gene_names.tsv")
genes <- genes$genes

# Define UI for cytoscope
ui <- fluidPage(
  
  # Application title
  titlePanel("Developing mouse brain"),
  
  # Sidebar with input
  sidebarLayout(
    sidebarPanel(width = 3,
                 selectInput("gene", "Gene", choices = genes, multiple = TRUE, selected = "Ttyh1"),
                 selectInput("mode", "How to plot expression", multiple = FALSE, choices = c("age", "pseudotime", "detected"), selected = "detected"),
                 numericInput("span", "LOESS smoothing", value = 1, min = 0.5, max = 1.5, step = 0.1),
                 numericInput("min_prop", "Minimum detection rate to show cell type", value = 0.1, min = 0, max = 1, step = 0.05),
                 selectInput("points", "Show points", choices = c(TRUE, FALSE), selected = TRUE),
                 selectInput("trend", "Show trend", choices = c(TRUE, FALSE), selected = TRUE),
                 selectInput("x_axis", "X-axis type", choices = c("even", "representative"), selected = "even"),
                 numericInput("jitter_width", "Width for point jitter", value = 0.3, min = 0, max = 0.5, step = 0.1)
              
                 
    ),
    
    # Output plots
    mainPanel(
      plotOutput("time", width = "10in", height = "4in"),
      plotOutput("time_legend"),
      downloadButton("download_time", "Download (png)")
    )
  )
)
