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
library(shinythemes)

options(shiny.maxRequestSize=6*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("yeti"),
  
  # Define title for the app and explain usage
  titlePanel(
    div(h2('Molecular mechanisms underlying plasticity in a thermally varying environment'), 
        h4("To use this application, first upload the sample information under the samples
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
                              min=0, max=100, value=75),
                 
                 sliderInput('slider_nonzero', 
                             'Select minimum non-zero samples:', 
                              min=0, max=80, value=20),
                 # slider to aggregate or dis-aggregate
                 #radioButtons('counts_ag', 
                 #            'Choose whether to disaggregate by male, female, larval:',
                 #              choices=c('Aggregate', 'Disaggregate'),
                 #              selected='Aggregate'),
                 
                 # plot button
                 actionButton('counts_button', 'Plot', icon= icon('chart-line'), width='100%'),
                 width=3
               ),
               
               # MAIN PANEL: GRAPH TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("Summary",
                            withSpinner(tableOutput("counts_summary"))),
                   tabPanel("Scatter Plots",
                            withSpinner(plotOutput("counts_scatter1")),
                            withSpinner(plotOutput("counts_scatter2"))),
                   tabPanel("Heatmap",
                            withSpinner(plotOutput("counts_heatmap"))),
                   tabPanel("PCA",
                            # sidebar with PCA input options
                            sidebarLayout(
                              sidebarPanel(
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
                                
                                selectInput('pca_color', 
                                            'Choose a color variable:',
                                            choices=c('Sex', 'Treatment', 'Timepoint'),
                                            selected='Treatment'),
                                
                                selectInput('pca_shape', 
                                            'Choose a symbol variable:',
                                            choices=c('Sex', 'Treatment', 'Timepoint'),
                                            selected='Sex'),
                                
                                # plot button
                                actionButton('pca_button', 'Update PCA', icon= icon('chart-line'), width='100%')
                              ),
                              mainPanel(
                                withSpinner(plotOutput("counts_pca")))
                              )
                            )
                   
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
                 
                 # options for comparison will depend on the radio button selection
                 uiOutput('de_comparison'),
                 
                 # slider
                 sliderInput('de_slider', 
                             'Select the magnitude of the p adjusted cutoff for each of the following analyses:', 
                             min=-80, max=0, value=-2),
                 
                 # plot button
                 actionButton('de_button', 'Get Dataset', icon= icon('chart-line'), width='100%'),
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
                                   'at the selected p-adjusted threshold. Note: adjusted p-values from LRT
                                   test are used.'),
                                
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
             sidebarLayout(
               # SIDEBAR: COUNTS DATA INPUT
               sidebarPanel(
                 # file input
                 h5(tags$b('Select a dataset to analyze. Use the sliders to specify a 
                           significance threshold and the minimum required genes per cluster.')),
                 # radio buttons
                 radioButtons('cluster_dataset', 
                              'Select a dataset for gene clustering:',
                              choices=c('Male', 'Female', 'Larvae'),
                              selected='Male'),
                 
                 # suggested pval and minc will change based on dataset chosen
                 uiOutput('cluster_pval_suggest'),
                 
                 # options for comparison will depend on the radio button selection
                 uiOutput('cluster_minc_suggest'),
                 
                 # plot button
                 actionButton('cluster_button', 'Plot', icon= icon('chart-line'), width='100%'),
                 width=3
               ),
               
               # MAIN PANEL: GRAPH TABS
               mainPanel(
                 headerPanel(""),
                 tabsetPanel(
                   tabPanel("Cluster Summary",
                            withSpinner(DT::dataTableOutput("cluster_info", width="100%"))),
                   tabPanel("Cluster Plots",
                            withSpinner(plotOutput("cluster_plots"))),
                   tabPanel("Timecourse Heatmap",
                            withSpinner(plotOutput("cluster_heatmap"))),
                   
                   tabPanel("Individual Gene",
                            # sidebar with gene input option
                            sidebarLayout(
                              sidebarPanel(
                                h5(tags$b('Track expression of a gene of interest for all replicates
                                          across timepoints.')),
                                
                                textInput('cluster_gene_input',
                                          'Enter an individual gene name:',
                                          width='100%',
                                          placeholder='ex. FBgn0003068'),
                                # plot button
                                actionButton('cluster_gene_button', 'Plot', icon= icon('chart-line'), width='100%')
                                
                              ),
                              mainPanel(
                                withSpinner(plotOutput("cluster_gene")))
                              )
                            )
                            
                 )
               ))
    ),

    ########## TAB 5 - UNDER DEVELOPMENT: GSEA ########## 
    #tabPanel("GSEA",
    #         # tableOutput('table'),
    #)
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  ### DEFINING UI OUTPUTS ###
  output$de_comparison <- renderUI({
    
    if(input$de_choice == 'Larvae'){
      elem <- radioButtons('de_comparison_choice', 
                           'Select a comparison:',
                           choices=c('Full vs. reduced model (LRT)', 
                                     'Fluctuating vs. Control at High1 (Wald)', 
                                     'Fluctuating vs. Control at Low1 (Wald)'),
                           selected='Full vs. reduced model (LRT)')
    } else{
      elem <- radioButtons('de_comparison_choice', 
                           'Select a comparison:',
                           choices=c('Full vs. reduced model (LRT)', 
                                     'Fluctuating vs. Control at High2 (Wald)', 
                                     'Fluctuating vs. Control at Low2 (Wald)',
                                     'Fluctuating vs. Control at High3 (Wald)', 
                                     'Fluctuating vs. Control at Low3 (Wald)'),
                           selected='Full vs. reduced model (LRT)')
    }
    elem
  })
  
  output$cluster_pval_suggest <- renderUI({
    if(input$cluster_dataset == 'Larvae'){
      elem <- sliderInput('cluster_pval', 
                          'Select the p-value threshold:', 
                          min=-55, max=0, value=-2)
      
    } else if(input$cluster_dataset == 'Male'){
      elem <- sliderInput('cluster_pval', 
                          'Select the p-value threshold:', 
                          min=-130, max=0, value=-5)
    } else if (input$cluster_dataset == 'Female'){
      elem <- sliderInput('cluster_pval', 
                  'Select the p-value threshold:', 
                  min=-70, max=0, value=-5)
    }
    elem
  })

  output$cluster_minc_suggest <- renderUI({
    if(input$cluster_dataset == 'Larvae'){
      elem <- sliderInput('cluster_minc', 
                          'Select the minimum number of genes per cluster:', 
                          min=0, max=50, value=10)
      
    } else if(input$cluster_dataset == 'Male'){
      elem <- sliderInput('cluster_minc', 
                          'Select the minimum number of genes per cluster:', 
                          min=0, max=100, value=30)
    } else if(input$cluster_dataset == 'Female'){
      elem <- sliderInput('cluster_minc', 
                  'Select the minimum number of genes per cluster:', 
                  min=0, max=100, value=30)
    }
    elem
  })
  
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
    # first, see if male, female, or larvae was chosen
    if (input$de_choice == 'Male'){
      id <- 'm'
    } else if (input$de_choice == 'Female'){
      id <- 'f'
    } else if (input$de_choice == 'Larvae'){
      id <- 'l'
    }
    
    # then, see what comparison was chosen
    if (input$de_comparison_choice == 'Full vs. reduced model (LRT)'){
      comp <- 'fluctvsctrl'
    } else if (input$de_comparison_choice == 'Fluctuating vs. Control at High1 (Wald)'){
      comp <- 'h1'
    } else if (input$de_comparison_choice == 'Fluctuating vs. Control at Low1 (Wald)'){
      comp <- 'l1'
    } else if (input$de_comparison_choice == 'Fluctuating vs. Control at High2 (Wald)'){
      comp <- 'h2'
    } else if (input$de_comparison_choice == 'Fluctuating vs. Control at Low2 (Wald)'){
      comp <- 'l2'
    } else if (input$de_comparison_choice == 'Fluctuating vs. Control at High3 (Wald)'){
      comp <- 'h3'
    } else if (input$de_comparison_choice ==  'Fluctuating vs. Control at Low3 (Wald)'){
      comp <- 'l3'
    }
      
    # get the user-selected DESeq2 results and change to proper format
    deseq_res <- readRDS(paste('objects/deseq_', id, '_', comp, sep=""))
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
  
  ### LOADING THE CLUSTERING DATA ###
  load_cluster_data <- reactive({
    # need to have the counts and sample data available
    req(input$sample_data)
    req(input$counts_data)
    
    # take a dependency on plot button
    input$cluster_button
    
    counts <- load_counts_data()
    meta <- load_sample_data()
    
    # load the selected dataset
    isolate({
      if (input$cluster_dataset == 'Male'){
        # divide samples up to prep for analysis
        male_titles <- meta[meta$'Sex.ch1' == 'male','title'] 
        # remove the male outlier
        male_titles <- male_titles[! male_titles %in% c('CAL2M')]
        # get individual counts
        male_counts <- counts[,male_titles]
        return(male_counts)
        
      } else if (input$cluster_dataset == 'Female'){
        female_titles <- meta[meta$'Sex.ch1' == 'female','title']
        female_counts <- counts[,female_titles]
        return(female_counts)
        
      } else if (input$cluster_dataset == 'Larvae'){
        larvae_titles <- meta[meta$'lifestage.ch1' == 'larvae','title']
        larvae_counts <- counts[,larvae_titles]
        return(larvae_counts)
      }
    })
    
  })
  
  ### LOADING THE CLUSTERING DE DATA ###
  # in this case we use fluct vs ctrl as default b/c we want the LRT comparison
  load_cluster_de_data <- reactive({
    if (input$cluster_dataset == 'Male'){
      deseq_res <- readRDS('objects/deseq_m_fluctvsctrl')
    } else if (input$cluster_dataset == 'Female'){
      deseq_res <- readRDS('objects/deseq_f_fluctvsctrl')
    } else if (input$cluster_dataset == 'Larvae'){
      deseq_res <- readRDS('objects/deseq_l_fluctvsctrl')
    }
    
    # get the user-selected DESeq2 results and change to proper format
    deseq_df <- as.data.frame(deseq_res)
    
    return(deseq_df)
    
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
      
      # CPM normalize
      counts <- normalize_by_cpm(as_tibble(counts, rownames='gene'))
      counts <-as.data.frame(counts)
      rownames(counts) <- counts$gene
      counts <- subset(counts, select = -c(gene))
      
      # make the table
      isolate({
      return(filter_counts(counts, input$slider_var/100, input$slider_nonzero))
      })
    }, caption="Summary of Genes Passing Filtering Conditions after CPM Normalization",
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
      
      # button for pca-specific functions
      input$pca_button
      
      counts <- load_counts_data()
      md <- load_sample_data()
      
      isolate({
        plot_pca(counts,
                 md,
                 input$counts_pcdim1,
                 input$counts_pcdim2,
                 input$pca_color,
                 input$pca_shape)
      })
      
    })
    
    ### CREATE EVENT REACTIVE FOR GETTING DESEQ2 DATASET ###
    get_deseq_selection <- eventReactive(input$de_button,{
      # require file inputs
      req(input$counts_data)
      req(input$sample_data)
      
      return(load_de_data())
      
    })
    
    ### MAKING DE DATA TABLE ###
    output$de_results <- DT::renderDataTable({
      # make it so input data is required
      req(input$counts_data)
      req(input$sample_data)
      
      # take a dependency on the plot button
      input$de_button
      
      # get the data
      de <- get_deseq_selection()
      
      isolate({
      
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
      
      # take dependency on volcano-specific button
      input$de_button
      
      data <- get_deseq_selection()
        
        volcano(data, 
                input$de_x_axis, 
                input$de_y_axis, 
                isolate({input$de_slider}),
                input$de_base,
                input$de_highlight)
      
    }, height=500)
    
    ### MAKING DE VENN DIAGRAM ###
    output$de_venn <- renderPlot({
      # make sure no error pops up before data uploaded
      req(input$counts_data)
      req(input$sample_data)
      
      # take dependency on input$plotbutton
      input$de_button
      
      # can't necessary compare all 3 for a given timepoint, so we'll just use full vs. reduced model
      data <- load_de_venn_data()
      
      # warn user instead of having error message pop up
      validate(need(length(data) >= 2, "Must select at least 2 datasets for venn diagram."))
      
      isolate({
      # get the de gene names (padj selected)
      gene_names <- lapply(data, de_genes, padj_cutoff=input$de_slider)

      ggVennDiagram(x=gene_names) +
        ggtitle(paste('Differentially Expressed Transcripts for padj<', input$de_slider),
                subtitle="Full vs. Reduced Model (LRT)") +
        scale_x_continuous(expand = expansion(mult = .2)) +
        theme(plot.title = element_text(face = "bold"))
      })
    })
    
    ### CREATE EVENT REACTIVE FOR DOING A NEW CLUSTERING THAT CAN BE USED ACROSS TABS ###
    make_new_clusters <- eventReactive(input$cluster_button,{
      # require file inputs
      req(input$counts_data)
      req(input$sample_data)
      
      # get data
      counts <- load_cluster_data()
      meta <- load_sample_data()
      deseq_results <- load_cluster_de_data()
      
      # do the clustering
      cluster_obj <- clustering(input$cluster_dataset, 
                                counts, 
                                meta, 
                                deseq_results, 
                                input$cluster_pval, 
                                input$cluster_minc)
      
    })

    ### GET SUMMARY OF CLUSTERING INFO ###
    output$cluster_info <- DT::renderDataTable(server=FALSE,{
      # make it so input data is required
      req(input$counts_data)
      req(input$sample_data)
      
      # get the cluster object
      cluster_obj <- make_new_clusters()
      
      # make data table
      DT::datatable(get_summary_data(cluster_obj),
                    extensions='Buttons',
                    filter='top',
                    options=list(
                    paging=TRUE,
                    searching=TRUE,
                    ordering=TRUE,
                    dom = 'Bfrtip',
                    columnDefs = list(
                      list(
                        targets = 0, className = "rownames"
                      )),
                    buttons = list(
                      list(extend = 'csv', 
                           filename='cluster_data',
                           exportOptions=list(columns = ":not(.rownames)"))),
                    pageLength=10))
      
    })
    
    ### GET CLUSTERING PLOTS ###
    output$cluster_plots <- renderPlot({
      # require file inputs
      req(input$counts_data)
      req(input$sample_data)
      
      # get the cluster object
      cluster_obj <- make_new_clusters()
      
      cluster_full_plot(cluster_obj)
      
    })
    
    ### GET TIMECOURSE HEATMAP ###
    output$cluster_heatmap <- renderPlot({
    # require file inputs
    req(input$counts_data)
    req(input$sample_data)
    
    # get the cluster object
    cluster_obj <- make_new_clusters()
    
    # make the heatmap
    timecourse_heatmap(cluster_obj)
    
    })
    
    ### CREATE EVENT REACTIVE FOR GETTING USER INPUT ###
    get_gene_name <- eventReactive(input$cluster_gene_button,{
      # require file inputs
      req(input$counts_data)
      req(input$sample_data)
      
      # get data
      return(input$cluster_gene_input)
      
    })
    
    ### GET INDIVIDUAL GENE PLOT ###
    output$cluster_gene <- renderPlot({
      # require file inputs
      req(input$counts_data)
      req(input$sample_data)
      
      # get the cluster object
      cluster_obj <- make_new_clusters()
      
      # get the user input and make sure it's valid first
      gene_name <- get_gene_name()
      validate(need(gene_name %in% cluster_obj$df$genes, 
                    "Gene name not included in clustering."))
      
      # get col data
      # take dependency on plot button
      input$cluster_gene_button
      
      isolate({
        
      counts <- load_cluster_data()
      meta <- load_sample_data()
      dataset <- input$cluster_dataset
        
      col_dat <- get_col_data(counts, meta, dataset)
      
      # make the plot
      cluster_gene_plot(cluster_obj, col_dat, gene_name)
        
      })
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
