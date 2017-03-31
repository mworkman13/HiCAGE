library(shiny)
library(HiCAGE)
library(tidyverse)

shinyUI(pageWithSidebar(
  headerPanel(list(em("HiCAGE:"), "Hi-C Annotation and Graphing Ensemble")),
  sidebarPanel(
    helpText(h4(list(em("Step 1"), "- Input Data Files:"))),
    fileInput('hicdata', 'Choose Hi-C data file'),
    fileInput('segmentdata', 'Choose Segmentation file'),
    fileInput('rnadata', 'Choose RNA-seq file'),
    helpText(h4(list(em("Step 2"), "- Select Data Columns:"))),
    textInput('hiccolumns',
              label = list(h5(strong("Enter Hi-C Column Numbers")),
                           h6("Order = Chrom1, Chrom1start, Chrom1end, Chrom2,
                              Chrom2start, Chrom2end, Score")),
              value = "1, 2, 3, 4, 5, 6, 8"),
    textInput('segmentcolumns',
              label = list(h5(strong("Enter Segmentation Column Numbers")),
                           h6("Order = Chromosome, Start, End, Mark, Score")),
              value = "1, 2, 3, 4, 5"),
    textInput('rnacolumns',
              label = list(h5(strong("Enter RNA Column Numbers")),
                           h6("Order = GeneID, FPKM")),
              value = c("1, 7")),
    helpText(h4(list(em("Step 3"), "- Select Species:"))),
    radioButtons("species", label = NULL,
                 choices = list("Human" = "hsapiens_gene_ensembl",
                                "Mouse" = "mmusculus_gene_ensembl"),
                 selected = "hsapiens_gene_ensembl"),
    helpText(h4(em("Step 4 - Go!"))),
    actionButton('goButton', "Run HiCAGE")
    ),
  mainPanel(
    tabsetPanel(id = "inTabset",
      tabPanel("Hi-C Data", dataTableOutput('contents')),
      tabPanel("Segmentation Data", dataTableOutput('contents2')),
      tabPanel("RNA-seq Data", dataTableOutput('contents3')),
      tabPanel("Final Data", dataTableOutput('finaldata'),
               downloadButton('downloadData', 'Download')),
      tabPanel("Circos Plot", plotOutput('circosplot')),
      tabPanel("Gene List", plotOutput('genelist'))
      )
    )
  ))
