#' Generates gene list for selected chromatin interactions
#'
#' @param datafile The GRanges object output from the overlap function
#' @param proximalmark The name of the mark in the GRanges object which you are
#' interested in obtaining nearest gene expression data from
#' @param distalmark The name of the mark in the GRanges object interacting
#' with the proximalmark
#' @param gene.symbol If TRUE, will ouput gene list using official gene symbol.
#' If FALSE, gene list output will use Ensembl ID
#' @param bio_mart BiomaRt used for pulling official gene name
#' @param martset The specific mart or dataset used for gene information
#' @param gbuild The genome version of mart to be used
#' @param geneOnto Run gene ontology enrichment using TopGO
#' @param expresscut Cutoff for expression value in gene list
#' @importFrom dplyr left_join bind_rows
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @import topGO
#' @import utils
#' @export
#' @return Returns a gene list sorted by gene expression levels of all genes
#' near the proximal mark that are interacting with the distal mark
#' @examples
#' hic_chr20 <- system.file("extdata", "hic_chr20.txt", package = "HiCAGE")
#' segment_chr20 <- system.file("extdata", "segment_chr20.bed",
#'     package = "HiCAGE")
#' rna_chr20 <- system.file("extdata", "rna_chr20.tsv", package = "HiCAGE")
#' overlapoutput <- overlap(hicfile = hic_chr20,
#'                          segmentfile = segment_chr20,
#'                          rnafile = rna_chr20)
#' gogenelist(datafile = overlapoutput, proximalmark = "PAR",
#'     distalmark = "EAR")

gogenelist <- function(datafile,
                       proximalmark,
                       distalmark,
                       gene.symbol = FALSE,
                       species = "human",
                       bio_mart = "ensembl",
                       martset = "hsapiens_gene_ensembl",
                       webhost = "www.ensembl.org",
                       geneOnto = TRUE,
                       expression_cutoff = 0) {
  finaldata <- as.data.frame(datafile)
  gene_list <- subset(finaldata, mark1==proximalmark & mark2==distalmark)
  gene_list$gene1 <- as.character(gene_list$gene1)
  gene_list$gene2 <- as.character(gene_list$gene2)
  gene_list <- subset(gene_list, select = c("gene1",
                                        "logFPKM1"))
  colnames(gene_list) <- c("gene", "logFPKM")
  gene_list2 <- subset(finaldata, mark1==distalmark & mark2==proximalmark)
  gene_list2$gene1 <- as.character(gene_list2$gene1)
  gene_list2$gene2 <- as.character(gene_list2$gene2)
  gene_list2 <- subset(gene_list2, select = c("gene2",
                                          "logFPKM2"))
  colnames(gene_list2) <- c("gene", "logFPKM")
  gene_list <- bind_rows(gene_list, gene_list2)
  gene_list <- gene_list[!duplicated(gene_list[c("gene",
                                           "logFPKM")]),]
  gene_list <- gene_list[order(-gene_list$logFPKM),]
  gene_list <- gene_list[!gene_list$logFPKM < expression_cutoff,]

  if (geneOnto == TRUE) {
    if (species == "human") {
      orgdb <- "org.Hs.eg.db"
    }
    if (species == "mouse") {
      orgdb <- "org.Mm.eg.db"
      }
    gene1 <- datafile[,c("gene1", "logFPKM1")]
    gene1 <- gene1[!gene1$logFPKM1 < expression_cutoff,]
    gene1 <- gene1[,1, drop = FALSE]
    colnames(gene1) [1] <- "gene"

    gene2 <- datafile[,c("gene2", "logFPKM2")]
    gene2 <- gene2[!gene2$logFPKM2 < expression_cutoff,]
    gene2 <- gene2[,1, drop = FALSE]
    colnames(gene2) [1] <- "gene"
    full_gene <- bind_rows(gene1, gene2)
    full_gene <- full_gene[!duplicated(full_gene),]
    full_gene <- as.factor(full_gene)

    inter_gene <- gene_list[,1]
    inter_gene <- as.factor(inter_gene)

    geneList <- factor(as.integer (full_gene %in% inter_gene))
    names(geneList) <- full_gene

    GOdata <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  description ="GO Analysis",
                  annot = annFUN.org,
                  mapping = orgdb,
                  ID = "Ensembl")
    resultFis <- runTest(GOdata,
                         algorithm = "classic",
                         statistic = "fisher")
    resultFisWeight <- runTest(GOdata,
                               algorithm = "weight",
                               statistic = "fisher")
    resultKS <- runTest(GOdata,
                        algorithm = "classic",
                        statistic = "ks")
    resultKS.elim <- runTest(GOdata,
                             algorithm = "elim",
                             statistic = "ks")
    allGO = usedGO(object = GOdata)
    allRes <- GenTable(GOdata,
                       classicFisher_p.value = resultFis,
                       weightedFisher_p.value = resultFisWeight,
                       classicKS_p.value = resultKS,
                       elimKS_p.value = resultKS.elim,
                       orderBy = "classicFisher_p.value",
                       topNodes = length(allGO))

  }
  else{}
  if (gene.symbol == TRUE) {
    ensembl = useMart(biomart = bio_mart,
                      dataset = martset,
                      host = webhost)
    gene <- getBM(attributes = c('ensembl_gene_id',
                                 'hgnc_symbol'),
                  mart = ensembl)

    gene_list <- left_join(gene_list, gene, by = c("gene" = "ensembl_gene_id"))
    gene_list <- subset(gene_list, select = c("gene",
                                              "hgnc_symbol",
                                              "logFPKM"))
  }
  else {}

  if (geneOnto == TRUE) {
    list('Gene_List' = gene_list, 'GO_Results' = allRes)

  }
  else{
    gene_list
  }
}
