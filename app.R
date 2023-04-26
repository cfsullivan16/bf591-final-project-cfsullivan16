## Caroline Sullivan
## BF591 Final Project
## Shiny App

# load functions from analysis script
source("analysis.R")

# load libraries
library(shiny)
library(DT)

options(shiny.maxRequestSize=6*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Define title for the app and explain usage
  titlePanel(
    div(h2('Molecular mechanisms underlying plasticity in a thermally varying environment'), 
        h4("To use this application..."))
  ),
  
  # Define 4-5 main tabs
  tabsetPanel(
    ########## TAB 1: SAMPLES ########## 
    tabPanel("Samples",
             sidebarLayout(
               # SIDEBAR: METADATA INPUT
               sidebarPanel(
                 # file input
                 fileInput('sample_data', 'Load sample info', accept=c('.csv')),
                 h5(tags$b('Select a sample variable to group by and a sample variable 
                    to count by in the visualizations tab.')),
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
                 actionButton('sample_button', 'Plot', icon= icon('chart-line'), width='100%')
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
    
    ########## TAB 2: COUNTS ########## 
    tabPanel("Counts",
             sidebarLayout(
               # SIDEBAR: COUNTS DATA INPUT
               sidebarPanel(
                 # file input
                 fileInput('counts_data', 'Load counts data', accept=c('.csv')),
                 h5(tags$b('Use the sliders to filter for genes above a certain variance
                     percentile and with at least a specified number of samples
                     with non-zero counts.')),
                 # radio buttons
                 # sliders
                 sliderInput('slider_var', 
                             'Select minimum percentile variance:', 
                              min=0, max=100, value=20),
                 
                 sliderInput('slider_nonzero', 
                             'Select minimum non-zero samples:', 
                              min=0, max=80, value=8),
                 
                 # plot button
                 actionButton('counts_button', 'Plot', icon= icon('chart-line'), width='100%')
               ),
               
               # MAIN PANEL: GRAPH TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("counts_summary")),
                   tabPanel("Scatter Plots",
                            plotOutput("counts_scatter1"),
                            plotOutput("counts_scatter2")),
                   tabPanel("Heatmap",
                            plotOutput("counts_heatmap")),
                   tabPanel("PCA",
                            plotOutput("counts_pca"))
               )
             ))
    ),
    ########## TAB 3: DE ########## 
    tabPanel("DE",
             # tableOutput('table'),
    ),
    ########## TAB 4: CLUSTERING ########## 
    tabPanel("Clustering",
             # tableOutput('table'),
    ),
    ########## TAB 5 IF TIME: GSEA ########## 
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
    file <- input$sample_data
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Must upload a csv file."))
    
    # get the metadata
    meta <- read.table(file$datapath, sep=',', header=TRUE, row.names = 1)
    
    return(meta)
  })
  
  ### LOADING THE COUNTS DATA ###
  load_counts_data <- reactive({
    # test that it is the correction data type, send error msg if not
    file <- input$counts_data
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Must upload a csv file."))
    
    # get the counts data
    counts <- read.csv(file$datapath)
    
    # turn gene_ids into rownames
    rownames(counts) <- counts$gene_id
    counts = subset(counts, select = -c(gene_id))
    
    return(counts)
  })
    
    
    ### MAKING SAMPLE SUMMARY TABLE ###
    output$sample_summary <- renderTable({
      # make sure no error pops up before data uploaded
      req(input$sample_data)
      
      # take a dependency on plot button
      input$sample_button
      
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
      input$sample_button
      
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
      input$sample_button
      
      # load the input data
      md <- load_sample_data()
      
      isolate({
        visualize_md(filter_md(md), 
                     input$group, 
                     input$metric)
      })
    }, width=500)
    
    ### MAKING COUNTS SUMMARY TABLE ###
    output$counts_summary <- renderTable({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      
      # take a dependency on plot button
      input$counts_button
      
      # load the data
      counts <- load_counts_data()
      
      # make the table
      isolate({
      return(filter_counts(counts, input$slider_var/100, input$slider_nonzero))
      })
    }, caption="Summary for Genes Passing Filtering Conditions",
       height=500) 
    
    ### MAKING COUNTS SCATTER PLOTS ###
    output$counts_scatter1 <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      
      # take dependency on button
      input$counts_button
      
      counts <- load_counts_data()
      
      isolate({
        scatter_plot1(prep_scatter_data(counts),
                      input$slider_var/100, 
                      input$slider_nonzero)
      })
      
    })
    
    output$counts_scatter2 <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      
      # take dependency on button
      input$counts_button
      
      counts <- load_counts_data()
      
      isolate({
        scatter_plot2(prep_scatter_data(counts),
                      input$slider_var/100, 
                      input$slider_nonzero)
      })
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
