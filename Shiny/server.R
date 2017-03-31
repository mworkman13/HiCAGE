options(shiny.maxRequestSize=100*1024^2)
shinyServer(function(input, output, session) {
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
                   hic.columns = as.numeric(unlist(strsplit(input$hiccolumns,","))),
                   segment.columns = as.numeric(unlist(strsplit(input$segmentcolumns,","))),
                   rna.columns = as.numeric(unlist(strsplit(input$rnacolumns,","))),
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
    circleplot(get.overlaps(),
               plot.subset = plotsub)
      })
    #Prepare Circos Plot for Download
    output$downloadPlot <- downloadHandler(
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
                   plot.subset = plotsub)
        dev.off()},
      contentType = 'image/png')
    
    
    })

})
