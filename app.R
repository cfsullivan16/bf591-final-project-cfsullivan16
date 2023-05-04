## Caroline Sullivan
## BF591 Final Project
## Shiny App

# load functions from analysis script
source("analysis.R")

# load libraries
library(shiny)
library(colourpicker)
library(DT)
library(shinycssloaders)

options(shiny.maxRequestSize=6*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Define title for the app and explain usage
  titlePanel(
    div(h2('Molecular mechanisms underlying plasticity in a thermally varying environment'), 
        h4("To use this application first upload the sample information under the samples
           tab and the counts information under the counts tab. Then, navigate through the
           tabs and subtabs to explore the results."))
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
                    to count by in the visualize tab.')),
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
                 actionButton('sample_button', 'Plot', icon= icon('chart-line'), width='100%'),
                 width=3
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
                 # slider to aggregate or dis-aggregate
                 #radioButtons('counts_ag', 
                 #            'Choose whether to disaggregate by male, female, larval:',
                 #              choices=c('Aggregate', 'Disaggregate'),
                 #              selected='Aggregate'),
                 numericInput(
                   'counts_pcdim1',
                   'Select Principal Component for the x-axis:',
                   1,
                   min = 1,
                   max = 80,
                   step = 1,
                 ),
                 
                 numericInput(
                   'counts_pcdim2',
                   'Select Principal Component for the y-axis:',
                   2,
                   min = 1,
                   max = 80,
                   step = 1,
                 ),
                 
                 # plot button
                 actionButton('counts_button', 'Plot', icon= icon('chart-line'), width='100%'),
                 width=3
               ),
               
               # MAIN PANEL: GRAPH TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("counts_summary")),
                   tabPanel("Scatter Plots",
                            withSpinner(plotOutput("counts_scatter1")),
                            withSpinner(plotOutput("counts_scatter2"))),
                   tabPanel("Heatmap",
                            withSpinner(plotOutput("counts_heatmap"))),
                   tabPanel("PCA",
                            withSpinner(plotOutput("counts_pca")))
               )
             ))
    ),
    ########## TAB 3: DE ########## 
    tabPanel("DE",
             sidebarLayout(
               # SIDEBAR: DE DATA
               sidebarPanel(
                 # instructions
                 h5(tags$b('Choose which dataset you would like to analyze:')),
                 
                 # buttons to choose dataset
                 radioButtons('de_choice', 
                             'Select DESeq2 results for male, female, or larvae datasets:',
                               choices=c('Male', 'Female', 'Larvae'),
                               selected='Male'),
                 
                 # slider
                 sliderInput('de_slider', 
                             'Select the magnitude of the p adjusted cutoff for each of the following analyses:', 
                             min=-80, max=0, value=-2),
                 
                 # plot button
                 actionButton('de_button', 'Plot', icon= icon('chart-line'), width='100%'),
                 width=3              
               ),
               
               # MAIN PANEL: DE TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("DE Results",
                            DT::dataTableOutput("de_results", width="100%")),
               
                   tabPanel("Vocano plot",
                            # sidebar within volcano plot area for volcano plot-specific input
                            sidebarLayout(
                              sidebarPanel(
                                # instructions
                                h5('A volcano plot can be generated with ',
                                   tags$b('log2 fold-change'),
                                   'on the x-axis and ',
                                   tags$b('p-adjusted'), 
                                   'on the y-axis.'),
                                
                                # selecting axes for volcano plot
                                selectInput('de_x_axis', 
                                            'Choose the column for the x-axis',
                                            choices=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'),
                                            selected='log2FoldChange'),
                                
                                selectInput('de_y_axis', 
                                            'Choose the column for the y-axis',
                                            choices=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'),
                                            selected='padj'),
                                
                                # color inputs
                                colourInput("de_base", 'base point color', "black"),
                                colourInput("de_highlight", 'highlight point color', "lightblue")
                                
                              ),
                            
                    mainPanel(
                            withSpinner(plotOutput("de_volcano"))
                            )
                        )
                    ),
                    
                   tabPanel("Venn diagrams",
                            # sidebar with venn diagram-specific options
                            sidebarLayout(
                              sidebarPanel(
                                
                                # instructions
                                h5('A venn diagram can be generated to find ',tags$b('overlapping significant genes'),
                                   'at the selected p-adjusted threshold.'),
                                
                                # checkboxes for venn
                                checkboxGroupInput('de_venn_choice',
                                                   'Select a combination of datasets:',
                                                   choices=c('Male', 'Female', 'Larvae'),
                                                   selected=c('Male', 'Female', 'Larvae'))
                                
                              ),
                              mainPanel(
                            withSpinner(plotOutput("de_venn")))
                              )
                   )
                 )
               ))
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
    counts <- subset(counts, select = -c(gene_id))
    
    return(counts)
  })
  
  ### LOADING THE DE DATA ###
  load_de_data <- reactive({
    if (input$de_choice == 'Male'){
      id <- 'm'
    } else if (input$de_choice == 'Female'){
      id <- 'f'
    } else if (input$de_choice == 'Larvae'){
      id <- 'l'
    }
  
    # get the user-selected DESeq2 results and change to proper format
    deseq_res <- readRDS(paste('objects/deseq_', id, '_fluctvsctrl', sep=""))
    deseq_df <- as.data.frame(deseq_res)
    
    return(deseq_df)
    
  })
  
  ### LOADING THE VENN DIAGRAM DE DATA ###
  load_de_venn_data <- reactive({
    
    data <- list('deseq_results')
    
    # collect all the necessary datasets
    for(choice in input$de_venn_choice) { 
      if (choice == 'Male'){
        id <- 'm'
      } else if (choice == 'Female'){
        id <- 'f'
      } else if (choice == 'Larvae'){
        id <- 'l'
      }
      
      # get the user-selected DESeq2 results and change to proper format
      deseq_res <- readRDS(paste('objects/deseq_', id, '_fluctvsctrl', sep=""))
      
      # add to list
      data <- append(data, deseq_res)
    }
    
    # isolate the deseq results and give list names 
    data_final <- data[2:length(data)]
    names(data_final) <- input$de_venn_choice
    
      
      # return list of deseq results
      return(data_final)
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
        scatter_plot1(counts,
                      input$slider_var/100, 
                      input$slider_nonzero)
      })
      
    }, width=725)
    
    output$counts_scatter2 <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      
      # take dependency on button
      input$counts_button
      
      counts <- load_counts_data()
      
      isolate({
        scatter_plot2(counts,
                      input$slider_var/100, 
                      input$slider_nonzero)
      })
      
    }, width=725)
    
    ### MAKING COUNTS HEATMAP ###
    output$counts_heatmap <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$sample_data)
      req(input$counts_data)
      
      # take dependency on button
      input$counts_button
      
      counts <- load_counts_data()
      md <- load_sample_data()
      
      isolate({
        plot_heatmap(counts,
                      md,
                      input$slider_var/100, 
                      input$slider_nonzero)
      })
      
    }, width=600, height=600)
    
    ### MAKING COUNTS PCA ###
    output$counts_pca <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$sample_data)
      req(input$counts_data)
      
      input$counts_button
      
      counts <- load_counts_data()
      md <- load_sample_data()
      
      isolate({
        plot_pca(counts,
                 md,
                 input$counts_pcdim1,
                 input$counts_pcdim2)
      })
      
    }, width=700, height=500)
    
    ### MAKING DE DATA TABLE ###
    output$de_results <- DT::renderDataTable({
      # make it so input data is required
      req(input$counts_data)
      req(input$sample_data)
      
      # take a dependency on the plot button
      input$de_button
      
      isolate({
      # load in data given selected dataset
      # need to isolate b/c radio button used in load_de_data()
      de <- load_de_data()
      
      de_filtered <- de %>%
        dplyr::filter(padj < 1 * 10^input$de_slider)
      
      # reformat some of the larger decimals
      de_filtered$pvalue <- formatC(de_filtered$pvalue, digits=6) 
      de_filtered$padj <- formatC(de_filtered$padj, digits=6) 
      
      # make genes a column so they are sortable
      Gene <- rownames(de_filtered)
      de_filtered <- cbind(Gene, data.frame(de_filtered, row.names=NULL))
      })
      
      # make data table
      DT::datatable(de_filtered, options=list(scrollX=TRUE))
      
    })
    
    ### MAKING DE VOLCANO PLOT ###
    output$de_volcano <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      req(input$sample_data)
      
      # take dependency on input$plotbutton
      input$de_button
      
      isolate({
        data <- load_de_data()
        
        volcano(data, 
                input$de_x_axis, 
                input$de_y_axis, 
                input$de_slider,
                input$de_base,
                input$de_highlight)
      })
      
    }, height=500)
    
    ### MAKING DE VENN DIAGRAM ###
    output$de_venn <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      req(input$sample_data)
      
      # take dependency on de diagram button
      input$de_button
      
      isolate({
        
      data <- load_de_venn_data()
      
      # warn user instead of having error message pop up
      validate(need(length(data) >= 2, "Must select at least 2 datasets for venn diagram."))
      
      # get the de gene names (padj selected)
      gene_names <- lapply(data, de_genes, padj_cutoff=input$de_slider)

      ggVennDiagram(x=gene_names) +
        ggtitle('Differentially Expressed Transcripts fluctuating vs. control') +
        scale_x_continuous(expand = expansion(mult = .2)) +
        theme(plot.title = element_text(face = "bold"))
      
      })
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
