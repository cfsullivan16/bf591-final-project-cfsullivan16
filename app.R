## Caroline Sullivan
## BF591 Final Project
## Shiny App

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Define title for the app and explain usage
  titlePanel(
    div(h2('Molecular mechanisms underlying plasticity in a thermally varying environment'), 
        h4("To use this application..."))
  ),
  
  # Define a main panel with 4-5 tabs
  mainPanel(
    tabsetPanel(
      # TAB 1: SAMPLES
      tabPanel("Samples",
               sidebarLayout(
                 # SIDEBAR: METADATA INPUT
                 sidebarPanel(
                   # file input
                   fileInput('file', 'Load sample info', accept='.csv')),
                 # MAIN PANEL: GRAPH TABS
                 mainPanel(
                   headerPanel(""),
                   tabsetPanel(
                     tabPanel("Summary"),
                     tabPanel("Table"))
                   )
               )),
      
      # TAB 2: COUNTS
      tabPanel("Counts",
               # tableOutput('table'),
      ),
      # TAB 3: DE
      tabPanel("DE",
               # tableOutput('table'),
      ),
      # TAB 4: CLUSTERING
      tabPanel("Clustering",
               # tableOutput('table'),
      ),
      # TAB 5 IF TIME: GSEA
      tabPanel("GSEA?",
               # tableOutput('table'),
      )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
