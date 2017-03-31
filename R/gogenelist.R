#' Generates gene list for selected chromatin interactions
#'
#' @param datafile The GRanges object output from the overlap function
#' @param proximalmark The name of the mark in the GRanges object which you are
#' interested in obtaining nearest gene expression data from
#' @param distalmark The name of the mark in the GRanges object interacting
#' with the proximalmark
#' @param gene.symbol If TRUE, will ouput gene list using official gene symbol.
#' If FALSE, gene list output will use Ensembl ID
#' @param mart BiomaRt used for pulling official gene name
#' @param martset The specific mart or dataset used for gene information
#' @param gbuild The genome version of mart to be used
#' @importFrom readr read_tsv read_csv
#' @importFrom dplyr left_join bind_rows
#' @importFrom biomaRt useEnsembl getBM
#' @export
#' @return Returns a gene list sorted by gene expression levels of all genes
#' near the proximal mark that are interacting with the distal mark
#' @examples
#' gogenelist(datafile = overlapoutput, proximalmark = "PAR", distalmark = "EAR")

gogenelist <- function(datafile,
                       proximalmark,
                       distalmark,
                       gene.symbol = FALSE,
                       mart = "ensembl",
                       martset = "hsapiens_gene_ensembl",
                       gbuild = 37) {
  finaldata <- as.data.frame(datafile)
  gorilla <- subset(finaldata, mark1==proximalmark & mark2==distalmark)
  gorilla$gene1 <- as.character(gorilla$gene1)
  gorilla$gene2 <- as.character(gorilla$gene2)
  gorilla <- subset(gorilla, select = c("gene1",
                                        "logFPKM1"))
  colnames(gorilla) <- c("gene", "logFPKM")
  gorilla2 <- subset(finaldata, mark1==distalmark & mark2==proximalmark)
  gorilla2$gene1 <- as.character(gorilla2$gene1)
  gorilla2$gene2 <- as.character(gorilla2$gene2)
  gorilla2 <- subset(gorilla2, select = c("gene2",
                                          "logFPKM2"))
  colnames(gorilla2) <- c("gene", "logFPKM")
  gorilla <- bind_rows(gorilla, gorilla2)
  gorilla <- gorilla[!duplicated(gorilla[c("gene",
                                           "logFPKM")]),]
  gorilla <- gorilla[order(-gorilla$logFPKM),]

  if (gene.symbol == TRUE) {
    ensembl = useEnsembl(biomart = mart,
                         dataset = martset,
                         GRCh = gbuild)
    gene <- getBM(attributes=c('ensembl_gene_id',
                               'hgnc_symbol',
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                  mart = ensembl)
    gorilla <- left_join(gorilla, gene, by = c("gene" = "ensembl_gene_id"))
    gorilla <- subset(gorilla, select = c("gene",
                                          "hgnc_symbol",
                                          "logFPKM"))
  }
  else {}

  gorilla
}
