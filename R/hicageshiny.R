#' Launch HiCAGE Shiny application
#'
#' @import shiny
#' @import shinydashboard
#' @importFrom readr read_tsv
#' @export
#' @return When run, function will launch the HiCAGE Shiny application - a GUI
#' for user-friendly implementation of HiCAGE
hicageshiny <- function() {
  options(shiny.maxRequestSize=100*1024^2)
  shinyApp(
    ui = shinyUI(dashboardPage(
      dashboardHeader(title = list(em("HiCAGE:"), "Hi-C Annotation and Graphics Ensemble"),
                      titleWidth = 480),
      dashboardSidebar(
        sidebarMenu(id = "inTabset",
                    menuItem("Data Input", tabName = "datainput", icon = icon("folder-open")),
                    menuItem("Hi-C Data", tabName = "hicdatatab", icon = icon("list")),
                    conditionalPanel("input.inTabset === 'hicdatatab'",
                                     textInput('hiccolumns',
                                               label = list(h5(strong("Enter Hi-C Column Numbers")),
                                                            h6("Order = Chrom1, Chrom1start, Chrom1end,
                                                               Chrom2, Chrom2start, Chrom2end, Score")),
                                               value = "1, 2, 3, 4, 5, 6, 8")),
                    menuItem("Segmentation Data", tabName = "segmentdatatab", icon = icon("list-alt")),
                    conditionalPanel("input.inTabset === 'segmentdatatab'",
                                     textInput('segmentcolumns',
                                               label = list(h5(strong("Enter Segmentation Column Numbers")),
                                                            h6("Order = Chromosome, Start, End,
                                                               Mark, Score")),
                                               value = "1, 2, 3, 4, 5")),
                    menuItem("RNA-seq Data", tabName = "rnaseqdatatab", icon = icon("list-ol")),
                    conditionalPanel("input.inTabset === 'rnaseqdatatab'",
                                     textInput('rnacolumns',
                                               label = list(h5(strong("Enter RNA Column Numbers")),
                                                            h6("Order = GeneID, FPKM")),
                                               value = c("1, 7"))),
                    menuItem("Final Data", tabName = "finaldatatab", icon = icon("th")),
                    menuItem("Circos Plot", tabName = "circosplottab", icon = icon("bullseye")),
                    menuItem("UpsetR Plot", tabName = "upsetplottab", icon = icon("bar-chart")),
                    menuItem("Gene List", tabName = "genelisttab", icon = icon("exchange")),
                    menuItem("GO Analysis", tabName = "goanalysistab", icon = icon("cogs")),
                    menuItem(actionButton('goButton', "Run HiCAGE"))
                                     )
                    ),
      dashboardBody(
        tags$head(
          tags$style(HTML('
                          .skin-blue .main-header .logo {
                          background-color: #3c8dbc;
                          }
                          .skin-blue .main-header .logo:hover {
                          background-color: #3c8dbc;
                          }
                          ')),
          tags$style(HTML(".dataTables_wrapper { overflow-x: scroll; }" )
          )),
        tabItems(
          tabItem(tabName = "datainput",
                  helpText(h4(list(em("Step 1"), "- Select Genome and Build:"))),
                  selectInput("genomes", label = NULL,
                              choices = list("Human GRCh37(hg19)" = 1,
                                             "Human GRCh38(hg38)" = 2,
                                             "Mouse NCBI Build 37(mm9)" = 3,
                                             "Mouse GRCm38(mm10)" = 4),
                              selected = 1),
                  helpText(h4(list(em("Step 2"), "- Input Data Files:"))),
                  helpText(h5("NOTE: Sample data can be found at http://www.github.com/mworkman13/HiCAGE/inst/extdata")),
                  fileInput('hicdata', 'Choose Hi-C data file: Use a tab-delimited *.txt file or tab separated value *.tsv file'),
                  fileInput('segmentdata', 'Choose Segmentation file: Use a *.bed file or tab-delimited *.txt file'),
                  fileInput('rnadata', 'Choose RNA-seq file: Use a tab-delimited *.txt file or tab separated value *.tsv file.
                            Currently, only Ensembl IDs are supported'),
                  helpText(h4(list(em("Step 3"), "- Make sure to select the appropriate 'Data Columns' and
                                   'Order' for the Hi-C, Segmentation, and RNA-seq files using the sidebar
                                   (default values may not be appropriate for your data structure)"))),
                  helpText(h4(list(em("Step 4"), "- Run", em("HiCAGE"))))
                  ),

          tabItem(tabName = "hicdatatab",
                  dataTableOutput('contents')),
          tabItem(tabName = "segmentdatatab",
                  dataTableOutput('contents2')),
          tabItem(tabName = "rnaseqdatatab",
                  dataTableOutput('contents3')),
          tabItem(tabName = "finaldatatab",
                  dataTableOutput('finaldata'),
                  downloadButton('downloadData', 'Download')),
          tabItem(tabName = "circosplottab",
                  plotOutput('circosplot'),
                  fluidRow(column(3, sliderInput("hicscale",
                                                 "Hi-C Score Range:",
                                                 min = 0,
                                                 max = 1000,
                                                 value = c(0, 100),
                                                 step = 10)),
                           column(4,
                                  sliderInput(
                                    "rnascale",
                                    "RNA Expression Range:",
                                    min = 0,
                                    max = 20,
                                    value = c(0, 5),
                                    step = 0.5)),
                           column(4,
                                  actionButton("reload", "Reload"))),
                  downloadButton('downloadCircos', 'Download')),
          tabItem(tabName = "upsetplottab",
                  plotOutput('upsetplot'),
                  tags$script("$(document).on('shiny:connected',
                              function(event) {
                              var myWidth = $(window).width();
                              Shiny.onInputChange('shiny_width',myWidth)});"),
                  tags$script("$(document).on('shiny:connected',
                              function(event) {
                              var myHeight = $(window).height();
                              Shiny.onInputChange('shiny_height',myHeight)});"),
                  downloadButton('downloadUpset', 'Download')),
          tabItem(tabName = "genelisttab",
                  fluidRow(column(3, uiOutput("Box1")),
                           column(4, uiOutput("Box2")),
                           column(5, h5(strong("Run GO Analysis with current set?"))),
                           column(5, actionButton("geneonto", "Run GO Analysis!"))),
                  fluidRow(column(3, sliderInput("genecut",
                                                 "Set Gene Expression Cutoff",
                                                 min = 0,
                                                 max = 20,
                                                 value = 0.1,
                                                 step = 0.1))),
                  dataTableOutput('genelist'),
                  downloadButton('downloadGene', 'Download')),
          tabItem(tabName = "goanalysistab",
                  dataTableOutput('godata'),
                  downloadButton('downloadgo', 'Download'))
                  )
                  )
        )),

    server = function(input, output, session) {  #Load Hi-C data file
      output$contents <- renderDataTable({
        inFile1 <- input$hicdata

        if (is.null(inFile1))
          return(NULL)

        read_tsv(inFile1$datapath,
                 comment = "#",
                 guess_max = 100000)
      })

      #Load Segmentation data file
      output$contents2 <- renderDataTable({
        inFile2 <- input$segmentdata

        if (is.null(inFile2))
          return(NULL)

        read_tsv(inFile2$datapath,
                 col_names = FALSE,
                 comment = "#",
                 skip = 1)
      })

      #Load RNA-seq data file
      output$contents3 <- renderDataTable({
        inFile3 <- input$rnadata

        if (is.null(inFile3))
          return(NULL)

        read_tsv(inFile3$datapath,
                 comment = "#")
      })

      #Render Final Data Frame
      observeEvent(input$goButton, {
        get.overlaps <- reactive({
          isolate({
            inFile1 <- input$hicdata
            inFile2 <- input$segmentdata
            inFile3 <- input$rnadata
          })

          ##Pick genome for biomaRt based on selected from dropbox
          if (input$genomes == 1) {
            genome <- "hsapiens_gene_ensembl"
            host <- "http://Feb2014.archive.ensembl.org"
          }
          else if (input$genomes == 2){
            genome <- "hsapiens_gene_ensembl"
            host <- "www.ensembl.org"
          }
          else if (input$genomes == 3){
            genome <- "mmusculus_gene_ensembl"
            host <- "http://Feb2014.archive.ensembl.org"
          }
          else {
            genome <- "mmusculus_gene_ensembl"
            host <- "www.ensembl.org"
          }

          ##Run the overlap function from HiCAGE package
          return(overlap(inFile1$datapath,
                         inFile2$datapath,
                         inFile3$datapath,
                         hic.columns =
                           as.numeric(unlist(strsplit(input$hiccolumns,","))),
                         segment.columns =
                           as.numeric(unlist(strsplit(
                             input$segmentcolumns,","))),
                         rna.columns =
                           as.numeric(unlist(strsplit(input$rnacolumns,","))),
                         martset = genome,
                         webhost = host))
        })

        #Show progress bar on GUI
        output$finaldata <- renderDataTable({
          withProgress(message = 'Processing',
                       detail = 'This may take a moment...', value = 0, {
                         result <- get.overlaps()
                         incProgress(1, detail = "Completed")
                       })
          return(result)
        })

        #Change to Final Data tab when Run button is pressed
        observeEvent(input$goButton, {
          updateTabsetPanel(session, "inTabset",
                            selected = "finaldatatab")
        })



        #Prepare Final table of all interaction for downloading to csv file
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("data-", Sys.time(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(get.overlaps(), file, row.names = FALSE)
          }
        )

        #Generate Circos Plot
        ##Set plot.subset to plot only Hi-C data if RNA-seq file is missing
        if (is.null(input$rnadata)) {
          plotsub <- "hicscore"
        }
        else {
          plotsub <- FALSE
        }


        circos.plot <- reactive({
          input$reload
          isolate({
            circleplot(get.overlaps(),
                       plot.subset = plotsub,
                       hic.range = c(input$hicscale),
                       rna.range = c(input$rnascale))
          })
        })
        output$circosplot <- renderPlot({
          withProgress(message = 'Processing',
                       detail = 'This may take a moment...', value = 0, {
                         result <- circos.plot()
                         incProgress(1, detail = "Completed")
                       })
          print(result)
        })
        ##Prepare Circos Plot for Download
        output$downloadCircos <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.time(), '.png', sep='')},
          content = function(file) {
            png(file,
                width = 9,
                height = 9,
                units = "in",
                pointsize = 20,
                bg = "white",
                res = 300)
            circleplot(get.overlaps(),
                       plot.subset = plotsub,
                       hic.range = c(input$hicscale),
                       rna.range = c(input$rnascale))
            dev.off()},
          contentType = 'image/png')

        #Generate UpSetR Plot
        output$upsetplot <- renderPlot({
          withProgress(message = 'Processing',
                       detail = 'This may take a moment...', value = 0, {
                         result <- plotup(get.overlaps())
                         incProgress(1, detail = "Completed")
                       })
          return(result)
        })

        ##Prepare UpSetR plot for download
        output$downloadUpset <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.time(), '.png', sep='')},
          content = function(file) {
            png(file,
                width = input$shiny_width,
                height = input$shiny_height,
                res = 150)
            plotup(get.overlaps())
            dev.off()},
          contentType = 'image/png')


        #Go Gene List
        if (is.null(input$rnadata)) {

        }
        else {
          output$Box1 = renderUI(
            selectInput("golist1",
                        label = "Get genes near:",
                        choices = c(unique(c(get.overlaps()$mark1,
                                             get.overlaps()$mark2)))))
          output$Box2 = renderUI(
            selectInput("golist2",
                        label = "Interacting with:",
                        choices = c(unique(c(get.overlaps()$mark1,
                                             get.overlaps()$mark2)))))

          output$genelist <- renderDataTable(
            gogenelist(get.overlaps(),
                       proximalmark = input$golist1,
                       distalmark = input$golist2,
                       geneOnto = FALSE,
                       expression_cutoff = input$genecut)
          )


          output$downloadGene <- downloadHandler(
            filename = function() {
              paste("genelist-",
                    Sys.time(), "-",
                    input$golist1, "-",
                    input$golist2, ".csv", sep="")
            },
            content = function(file) {
              write.csv(get.genelist(), file, row.names = FALSE)
            }
          )

          #GO Analysis
          ##Pick genome for biomaRt based on selected from dropbox
          if (input$genomes == 1) {
            gospecies <- "human"
          }
          else if (input$genomes == 2){
            gospecies <- "human"
          }
          else if (input$genomes == 3){
            gospecies <- "mouse"
          }
          else {
            gospecies <- "mouse"
          }
          go.analysis <- eventReactive(input$geneonto, {
            gogenelist(get.overlaps(),
                       proximalmark = input$golist1,
                       distalmark = input$golist2,
                       species = gospecies,
                       geneOnto = TRUE,
                       expression_cutoff = input$genecut)

          })
          output$godata <- renderDataTable({
            withProgress(message = 'Processing',
                         detail = 'This may take a moment...', value = 0, {
                           result <- go.analysis()$GO_Results
                           incProgress(1, detail = "Completed")
                         })
            return(result)
          })

          #Change to GO Analysis tab when Run button is pressed
          observeEvent(input$geneonto, {
            updateTabsetPanel(session, "inTabset",
                              selected = "goanalysistab")
          })

          output$downloadgo <- downloadHandler(
            filename = function() {
              paste("goanalysis-",
                    Sys.time(), "-",
                    input$golist1, "-",
                    input$golist2, ".csv", sep="")
            },
            content = function(file) {
              write.csv(go.analysis()$GO_Results, file, row.names = FALSE)
            }
          )

        }
      })
    })
}
