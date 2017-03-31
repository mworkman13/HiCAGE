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
      inFile3 <- input$rnadata
      inFile2 <- input$segmentdata
      })

    return(overlap(inFile1$datapath,
                   inFile2$datapath,
                   inFile3$datapath,
                   hic.columns = as.numeric(unlist(strsplit(input$hiccolumns,","))),
                   segment.columns = as.numeric(unlist(strsplit(input$segmentcolumns,","))),
                   rna.columns = as.numeric(unlist(strsplit(input$rnacolumns,","))),
                   martset = input$species))
    })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(get.overlaps(), file, row.names = FALSE)
    }
  )
    output$finaldata <- renderDataTable({
      withProgress(message = 'Processing',
                   detail = 'This may take a moment...', value = 0, {
                     result <- get.overlaps()
                     incProgress(1, detail = "Completed")
                   })
      return(result)
    })

    updateTabsetPanel(session, "inTabset",
                      selected = "Final Data")
    #Generate Circos Plot
    output$circosplot <- renderPlot({
    circleplot(get.overlaps())
      })
    })

})
