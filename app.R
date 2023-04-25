## Caroline Sullivan
## BF591 Final Project
## Shiny App

# load functions from analysis script
source("analysis.R")

# load libraries
library(shiny)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Define title for the app and explain usage
  titlePanel(
    div(h2('Molecular mechanisms underlying plasticity in a thermally varying environment'), 
        h4("To use this application..."))
  ),
  
  # Define 4-5 main tabs
  tabsetPanel(
    # TAB 1: SAMPLES
    tabPanel("Samples",
             sidebarLayout(
               # SIDEBAR: METADATA INPUT
               sidebarPanel(
                 # file input
                 fileInput('sample_data', 'Load sample info', accept='.csv'),
                 # radio buttons
                 radioButtons('group', 
                             'Choose a category to group by',
                              choices=c('Lifestage', 'Sex', 'Timepoint', 'Treatment'),
                              selected='Lifestage'),
               
                 radioButtons('metric', 
                             'Choose a category to plot',
                              choices=c('Lifestage', 'Sex', 'Timepoint', 'Treatment'),
                              selected='Timepoint'),
                 # plot button
                 actionButton('plot_samples_button', 'Plot', icon= icon('chart-line'), width='100%')
                 ),
               
               # MAIN PANEL: GRAPH TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("sample_summary")),
                   tabPanel("Table",
                            DT::dataTableOutput("sample_info", width="100%")),
                   tabPanel("Visualize",
                            plotOutput("sample_plot")))
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


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ### LOADING THE SAMPLE DATA ###
  load_sample_data <- reactive({
    # test that it is the correction data type, send error msg if not
    
    
    # get the metadata
    meta <- read.table(input$sample_data$datapath, sep=',', header=TRUE, row.names = 1)
    
    return(meta)
  })
    
    
    ### MAKING SAMPLE SUMMARY TABLE ###
    output$sample_summary <- renderTable({
      # make sure no error pops up before data uploaded
      req(input$sample_data)
      
      # take a dependency on plot button
      input$plot_samples_button
      
      # load the data
      meta <- load_sample_data()
      
      # make the table
      return(summarize_md(filter_md(meta)))
    }, width="100%") 
    
    
    ### MAKING SAMPLE DATA TABLE ###
    output$sample_info <- DT::renderDataTable({
      # make it so input data is required
      req(input$sample_data)
      
      # take a dependency on the plot button
      input$plot_samples_button
      
      # load in data
      md <- load_sample_data()
      
      # make data table
      DT::datatable(filter_md(md),  options=list(scrollX=TRUE, pageLength=5))
      
    })
    
    ### MAKING SAMPLE VISUALIZATIONS ###
    output$sample_plot <- renderPlot({
      # require input data
      req(input$sample_data)
      
      # take dependency on plot button
      input$plot_samples_button
      
      # load the input data
      md <- load_sample_data()
      
      isolate({
        visualize_md(filter_md(md), 
                     input$group, 
                     input$metric)
      })
    }, width=500)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
