#' Launch HiCAGE Shiny application
#'
#' @import shiny
#' @importFrom readr read_tsv
#' @export
#' @return When run, function will launch the HiCAGE Shiny application - a GUI
#' for user-friendly implementation of HiCAGE
hicageshiny <- function() {
  options(shiny.maxRequestSize=100*1024^2)
  shinyApp(
    ui = shinyUI(pageWithSidebar(
      headerPanel(list(em("HiCAGE:"), "Hi-C Annotation and Graphing Ensemble")),
      sidebarPanel(
        helpText(h4(list(em("Step 1"), "- Select Genome and Build:"))),
        selectInput("genomes", label = NULL,
                    choices = list("Human GRCh37(hg19)" = 1,
                                   "Human GRCh38(hg38)" = 2,
                                   "Mouse NCBI Build 37(mm9)" = 3,
                                   "Mouse GRCm38(mm10)" = 4), selected = 1),
        helpText(h4(list(em("Step 2"), "- Input Data Files:"))),
        fileInput('hicdata', 'Choose Hi-C data file'),
        fileInput('segmentdata', 'Choose Segmentation file'),
        fileInput('rnadata', 'Choose RNA-seq file'),
        helpText(h4(list(em("Step 3"), "- Select Data Columns:"))),
        textInput('hiccolumns',
                  label = list(h5(strong("Enter Hi-C Column Numbers")),
                               h6("Order = Chrom1, Chrom1start, Chrom1end,
                                  Chrom2, Chrom2start, Chrom2end, Score")),
                  value = "1, 2, 3, 4, 5, 6, 8"),
        textInput('segmentcolumns',
                  label = list(h5(strong("Enter Segmentation Column Numbers")),
                               h6("Order = Chromosome, Start, End,
                                  Mark, Score")),
                  value = "1, 2, 3, 4, 5"),
        textInput('rnacolumns',
                  label = list(h5(strong("Enter RNA Column Numbers")),
                               h6("Order = GeneID, FPKM")),
                  value = c("1, 7")),
        helpText(h4(em("Step 4 - Go!"))),
        actionButton('goButton', "Run HiCAGE")
                  ),
      mainPanel(
        tabsetPanel(id = "inTabset",
                    tabPanel("Hi-C Data",
                             dataTableOutput('contents')),
                    tabPanel("Segmentation Data",
                             dataTableOutput('contents2')),
                    tabPanel("RNA-seq Data",
                             dataTableOutput('contents3')),
                    tabPanel("Final Data",
                             dataTableOutput('finaldata'),
                             downloadButton('downloadData', 'Download')),
                    tabPanel("Circos Plot",
                             plotOutput('circosplot'),
                             fluidRow(column(3, sliderInput("hicscale",
                                                            "Hi-C Score Range:",
                                                            min = 0,
                                                            max = 1000,
                                                            value = c(0, 100))),
                                      column(4,
                                             sliderInput(
                                               "rnascale",
                                               "RNA Expression Range:",
                                               min = 0,
                                               max = 20,
                                               value = c(0, 5))),
                                      column(4,
                                             actionButton("reload", "Reload"))),
                             downloadButton('downloadCircos', 'Download')),
                    tabPanel("UpSetR Plot",
                             plotOutput('upsetplot'),
                             fluidRow(column(3,
                                             numericInput("numwidth",
                                                          label = h6("Download Width (inches)"),
                                                          value = 1)),
                                      (column(4,
                                              numericInput("numheight",
                                                           label = h6("Download Height (inches)"),
                                                           value = 1)))),
                             downloadButton(column(4,'downloadUpset', 'Download'))),
                    tabPanel("Gene List",
                             fluidRow(column(3, uiOutput("Box1")),
                                      column(4, uiOutput("Box2"))),
                             dataTableOutput('genelist'),
                             downloadButton('downloadGene', 'Download'))

        )
      )
      )),

    server = function(input, output, session) {
      #Load Hi-C data file
      output$contents <- renderDataTable({
        inFile1 <- input$hicdata

        if (is.null(inFile1))
          return(NULL)

        read_tsv(inFile1$datapath,
                 col_names = FALSE,
                 comment = "#",
                 skip = 1,
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
                 col_names = FALSE,
                 comment = "#",
                 skip = 1)
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
            build <- 37
          }
          else if (input$genomes == 2){
            genome <- "hsapiens_gene_ensembl"
            build <- NULL
          }
          else if (input$genomes == 3){
            genome <- "mmusculus_gene_ensembl"
            build <- 37
          }
          else {
            genome <- "mmusculus_gene_ensembl"
            build <- NULL
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
                         gbuild = build))
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
        updateTabsetPanel(session, "inTabset",
                          selected = "Final Data")

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

        output$circosplot <- renderPlot({
          input$reload
          isolate({
          circleplot(get.overlaps(),
                     plot.subset = plotsub,
                     hic.range = c(input$hicscale),
                     rna.range = c(input$rnascale))
        })
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
          plotup(get.overlaps())
        })
        ##Prepare UpSetR plot for download
        output$downloadUpset <- downloadHandler(
          filename = function() {
            paste("plot-", Sys.time(), '.png', sep='')},
          content = function(file) {
            png(file,
                width = input$numwidth,
                height = input$numheight,
                units = "in",
                pointsize = 20,
                bg = "white",
                res = 300)
            plotup(get.overlaps())
            dev.off()},
          contentType = 'image/png')


        #Go Gene List
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

        get.genelist <- reactive({
          return(gogenelist(get.overlaps(),
                            proximalmark = input$golist1,
                            distalmark = input$golist2))
        })

        output$genelist <- renderDataTable(
          return(get.genelist()))

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

      })

    })
}
