

gogenelist <- function(datafile,
                       proximalmark,
                       distalmark,
                       gene.symbol = FALSE) {
  finaldata <- as.data.frame(datafile)
  gorilla <- subset(finaldata, mark1==proximalmark & mark2==distalmark)
  gorilla$gene1 <- as.character(gorilla$gene1)
  gorilla$gene2 <- as.character(gorilla$gene2)
  gorilla <- subset(gorilla, select = c("gene1",
                                        "FPKM1"))
  colnames(gorilla) <- c("gene", "FPKM")
  gorilla2 <- subset(finaldata, mark1==distalmark & mark2==proximalmark)
  gorilla2$gene1 <- as.character(gorilla2$gene1)
  gorilla2$gene2 <- as.character(gorilla2$gene2)
  gorilla2 <- subset(gorilla2, select = c("gene2",
                                          "FPKM2"))
  colnames(gorilla2) <- c("gene", "FPKM")
  gorilla <- bind_rows(gorilla, gorilla2)
  gorilla <- gorilla[!duplicated(gorilla[c("gene",
                                           "FPKM")]),]
  gorilla <- gorilla[order(-gorilla$FPKM),]

  if (gene.symbol == TRUE) {
    ensembl = useEnsembl(biomart="ensembl",
                         dataset="hsapiens_gene_ensembl",
                         GRCh=37)
    gene <- getBM(attributes=c('ensembl_gene_id',
                               'hgnc_symbol',
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                  mart = ensembl)
    gorilla <- left_join(gorilla, gene, by = c("gene" = "ensembl_gene_id"))
    gorilla <- subset(gorilla, select = c("gene",
                                          "hgnc_symbol",
                                          "FPKM"))
  }
  gorilla
}
