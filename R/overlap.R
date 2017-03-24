#' Combine chromosome conformation capture data with segmentation and RNA data
#'
#' @param hicfile The chromosome conformation capture datafile
#' @param segmentfile The segmentation datafile from StateHub
#' @param rnafile The RNA gene expression datafile containing Ensembl IDs and
#' FPKM or TPM normalized expression values
#' @param manualselect If TRUE, allows for manual selection of chromosome
#' conformation capture columns to be used for analysis
#' @param mart BiomaRt used for pulling gene name and location data
#' @param martset The specific mart or dataset used for genome information
#' @param gbuild The genome version of mart to be used
#' @param hic.columns Columns from the chromsome conformation capture data
#' that contains, in order, 'left chromosome', 'left start', 'left end', 'right
#' chromosome', 'right start', 'right end', 'interaction score'
#' @importFrom readr read_tsv read_csv
#' @importFrom tidyr separate
#' @importFrom dplyr mutate left_join full_join group_by
#' @importFrom GenomicRanges GRanges nearest findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom magrittr %>%
#' @export
#' @return Returns a GRanges with each chromsome conformation capture segment
#' annotated with a prioritized segmentation mark, the nearest gene to the
#' segmentation mark and the gene's FPKM expression data
#' @examples
#' overlap("hicfile.txt", "segmentfile.bed", "rnafile.tsv")
overlap <- function(hicfile,
                    segmentfile,
                    rnafile,
                    mart = "ensembl",
                    martset = "hsapiens_gene_ensembl",
                    gbuild = 37,
                    hic.columns = c(1:6,8),
                    segment.columns = c(1:4),
                    rna.columns = c(1,8)) {
  # Load biomart for assigning nearest gene and subset in GRanges
  ensembl = useEnsembl(biomart = mart,
                       dataset = martset,
                       GRCh = gbuild)
  gene <- getBM(attributes = c('ensembl_gene_id',
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                mart = ensembl)

  gene$chromosome_name <- gsub("^chr", "", gene$chromosome_name)
  ensGene <- GRanges(seqnames = gene$chromosome_name,
                     ranges = IRanges(start = as.integer(gene$start_position),
                                      end = as.integer(gene$end_position)))
  # Manual select allows user to manual choose columns from Hi-C datafile
    HiCdata <- read.table(file = hicfile,
                          comment.char = "#",
                          skip = 1)
    HiCdata <- subset(HiCdata, select = hic.columns)
    colnames(HiCdata) <- c("region1chrom",
                           "region1start",
                           "region1end",
                           "region2chrom",
                           "region2start",
                           "region2end",
                           "score")
  #Remove duplicates from Hi-C data and subset into GRanges
  HiCdata <- HiCdata[!duplicated(HiCdata[c("region1start",
                                           "region1end",
                                           "region2start",
                                           "region2end")]),]

  HiCdata$region1chrom <- gsub("^chr", "", HiCdata$region1chrom)
  HiCdata$region2chrom <- gsub("^chr", "", HiCdata$region2chrom)

  chiaregion1 <- subset(HiCdata, select = c(1:3))
  chiaregion1 <- GRanges(seqnames = chiaregion1$region1chrom,
                         ranges = IRanges(start = as.integer(chiaregion1$region1start),
                                          end = as.integer(chiaregion1$region1end)))

  chiaregion2 <- subset(HiCdata, select = c(4:6))
  chiaregion2 <- GRanges(seqnames = chiaregion2$region2chrom,
                         ranges = IRanges(start = as.integer(chiaregion2$region2start),
                                          end = as.integer(chiaregion2$region2end)))

  #Parse Segmentation data and subset in GenomicRanges
  #Segmentation files from StateHub Default Model
  #http://statehub.org/modeltracks/default_model/
  segmentation <- read_tsv(file = segmentfile,
                           col_names = FALSE,
                           comment = "#",
                           skip = 1)
  segmentation <- subset(segmentation, select = segment.columns)
  colnames(segmentation) <- c("chromosome",
                              "segstart",
                              "segend",
                              "state")

  segmentation$chromosome <- gsub("^chr", "", segmentation$chromosome)
  epigenetic <- GRanges(seqnames = segmentation$chromosome,
                        ranges = IRanges(start = as.integer(segmentation$segstart),
                                         end = as.integer(segmentation$segend)),
                        mark = as.character(segmentation$state))
  #Parse RNA seq file
  RNAseq <- read_tsv(rnafile) %>%
    separate("gene_id", into = c("gene_id", "extraint"), sep = "\\.")
  RNAseq <- subset(RNAseq, select = rna.columns)

  #Find nearest gene to each segmentation and Hi-C location; Add FPKM
  epineargene <- data.frame(nearest(epigenetic, ensGene, ignore.strand=TRUE))
  colnames(epineargene) <- "ensembl"
  epineargene <- data.frame(gene[epineargene$ensembl, "ensembl_gene_id"])
  colnames(epineargene) <- "ensembl"
  epineargene$ensembl <- as.character(epineargene$ensembl)
  epineargene <- left_join(epineargene, RNAseq, by = c("ensembl" = "gene_id"))
  segmentation <- cbind(segmentation, epineargene)

  HiCneargene1 <- data.frame(nearest(chiaregion1, ensGene, ignore.strand=TRUE))
  colnames(HiCneargene1) <- "ensembl"
  HiCneargene1 <- data.frame(gene[HiCneargene1$ensembl, "ensembl_gene_id"])
  colnames(HiCneargene1) <- "ensembl"
  HiCneargene1$ensembl <- as.character(HiCneargene1$ensembl)
  HiCneargene1 <- left_join(HiCneargene1, RNAseq, by = c("ensembl" = "gene_id"))

  HiCneargene2 <- data.frame(nearest(chiaregion2, ensGene, ignore.strand=TRUE))
  colnames(HiCneargene2) <- "ensembl"
  HiCneargene2 <- data.frame(gene[HiCneargene2$ensembl, "ensembl_gene_id"])
  colnames(HiCneargene2) <- "ensembl"
  HiCneargene2$ensembl <- as.character(HiCneargene2$ensembl)
  HiCneargene2 <- left_join(HiCneargene2, RNAseq, by = c("ensembl" = "gene_id"))

  #Intersect ChIA-PET data with Epigenetic Segmentation Data

  overlap1 <- data.frame(findOverlaps(chiaregion1, epigenetic))

  overlap2 <- data.frame(findOverlaps(chiaregion2, epigenetic))

  region1data <- data.frame(overlap1[, "queryHits"],
                            HiCdata[overlap1$queryHits, "region1chrom"],
                            HiCdata[overlap1$queryHits, "region1start"],
                            HiCdata[overlap1$queryHits, "region1end"],
                            segmentation[overlap1$subjectHits,"state"],
                            segmentation[overlap1$subjectHits,"ensembl"],
                            segmentation[overlap1$subjectHits, "FPKM"])
  colnames(region1data) <- c("overlap1....queryHits..",
                             "region1chrom",
                             "region1start",
                             "region1end",
                             "state",
                             "ensembl1",
                             "FPKM1")

  region2data <- data.frame(overlap2[, "queryHits"],
                            HiCdata[overlap2$queryHits, "region2chrom"],
                            HiCdata[overlap2$queryHits, "region2start"],
                            HiCdata[overlap2$queryHits, "region2end"],
                            segmentation[overlap2$subjectHits,"state"],
                            segmentation[overlap2$subjectHits,"ensembl"],
                            segmentation[overlap2$subjectHits, "FPKM"])
  colnames(region2data) <- c("overlap2....queryHits..",
                             "region2chrom",
                             "region2start",
                             "region2end",
                             "state",
                             "ensembl2",
                             "FPKM2")

  #Prioritize epigenetic marks and remove multiple entries for a given segment
  region1data$state <- ordered(region1data$state,
                               levels = c("TRS",
                                          "SCR",
                                          "HET",
                                          "RPS",
                                          "ERC",
                                          "ER",
                                          "PRC",
                                          "PR",
                                          "EPRC",
                                          "EPR",
                                          "PPRC",
                                          "PPR",
                                          "CTCFC",
                                          "CTCF",
                                          "ARC",
                                          "AR",
                                          "EARC",
                                          "EAR",
                                          "PARC",
                                          "PAR"))
  region1data <- region1data %>%
    group_by(overlap1....queryHits..) %>%
    mutate(mark1=state[which.max(state)])

  region1data <- region1data[!duplicated
                             (region1data[c("overlap1....queryHits..")]),]

  results1 <- data.frame(HiCdata[, "region1chrom"],
                         HiCdata[, "region1start"],
                         HiCdata[, "region1end"])
  results1 <- results1 %>%
    mutate(ID=rownames(results1))
  results1$ID = as.integer(results1$ID)
  colnames(results1) <- c("region1chrom", "region1start", "region1end", "ID")
  results1 <- results1[c('ID', 'region1chrom', 'region1start', 'region1end')]

  region2data$state <- ordered(region2data$state,
                               levels = c("TRS",
                                          "SCR",
                                          "HET",
                                          "RPS",
                                          "ERC",
                                          "ER",
                                          "PRC",
                                          "PR",
                                          "EPRC",
                                          "EPR",
                                          "PPRC",
                                          "PPR",
                                          "CTCFC",
                                          "CTCF",
                                          "ARC",
                                          "AR",
                                          "EARC",
                                          "EAR",
                                          "PARC",
                                          "PAR"))
  region2data <- region2data %>%
    group_by(overlap2....queryHits..) %>%
    mutate(mark2=state[which.max(state)])

  region2data <- region2data[!duplicated
                             (region2data[c("overlap2....queryHits..")]),]

  results2 <- data.frame(HiCdata[, "region2chrom"],
                         HiCdata[, "region2start"],
                         HiCdata[, "region2end"])
  results2 <- results2 %>%
    mutate(ID=rownames(results2))
  results2$ID = as.integer(results2$ID)
  colnames(results2) <- c("region2chrom", "region2start", "region2end", "ID")
  results2 <- results2[c('ID', 'region2chrom', 'region2start', 'region2end')]

  #Combine HiCdata data with prioritized epigenetic mark
  results1 <- left_join(results1,
                        region1data[, c("overlap1....queryHits..",
                                        "mark1",
                                        "ensembl1",
                                        "FPKM1")],
                        by = c("ID" = "overlap1....queryHits.."))


  results2 <- left_join(results2,
                        region2data[, c("overlap2....queryHits..",
                                        "mark2",
                                        "ensembl2",
                                        "FPKM2")],
                        by = c("ID" = "overlap2....queryHits.."))

  #Replace NA genes and FPKM with nearest gene to HiC range
  results1$ensembl1 <- as.character(results1$ensembl1)
  results1$ensembl1[is.na(results1$ensembl1)] <-
    HiCneargene1$ensembl[is.na(results1$ensembl1)]

  results1$FPKM1[is.na(results1$FPKM1)] <-
    HiCneargene1$FPKM[is.na(results1$FPKM1)]

  results2$ensembl2 <- as.character(results2$ensembl2)
  results2$ensembl2[is.na(results2$ensembl2)] <-
    HiCneargene2$ensembl[is.na(results2$ensembl2)]

  results2$FPKM2[is.na(results2$FPKM2)] <-
    HiCneargene2$FPKM[is.na(results2$FPKM2)]

  #Combine results for HiC region 1 and region 2
  final <- data.frame(results1[,"ID"],
                      results1[, "region1chrom"],
                      results1[, "region1start"],
                      results1[, "region1end"],
                      results1[, "mark1"],
                      results1[, "ensembl1"],
                      results1[, "FPKM1"],
                      results2[, "region2chrom"],
                      results2[, "region2start"],
                      results2[, "region2end"],
                      results2[, "mark2"],
                      results2[, "ensembl2"],
                      results2[, "FPKM2"],
                      HiCdata[, "score"])

  final <- setNames(final, c("ID",
                             "region1chrom",
                             "region1start",
                             "region1end",
                             "mark1",
                             "gene1",
                             "FPKM1",
                             "region2chrom",
                             "region2start",
                             "region2end",
                             "mark2",
                             "gene2",
                             "FPKM2",
                             "score"))

  #Rename epigenetic marks in final data frame
  final$mark1 <- as.character.factor (final$mark1)
  final$mark2 <- as.character.factor (final$mark2)

  final[is.na(final)] <- "None"
  final$mark1[final$mark1=="AR"] <- "AR"
  final$mark1[final$mark1=="ARC"] <- "AR"
  final$mark1[final$mark1=="CTCF"] <- "CTCF"
  final$mark1[final$mark1=="CTCFC"] <- "CTCF"
  final$mark1[final$mark1=="EAR"] <- "EAR"
  final$mark1[final$mark1=="EARC"] <- "EAR"
  final$mark1[final$mark1=="EPR"] <- "EPR"
  final$mark1[final$mark1=="EPRC"] <- "EPR"
  final$mark1[final$mark1=="ER"] <- "ER"
  final$mark1[final$mark1=="ERC"] <- "ER"
  final$mark1[final$mark1=="HET"] <- "HET"
  final$mark1[final$mark1=="PAR"] <- "PAR"
  final$mark1[final$mark1=="PARC"] <- "PAR"
  final$mark1[final$mark1=="PPR"] <- "PPR"
  final$mark1[final$mark1=="PPRC"] <- "PPR"
  final$mark1[final$mark1=="PR"] <- "PR"
  final$mark1[final$mark1=="PRC"] <- "PR"
  final$mark1[final$mark1=="RPS"] <- "RPS"
  final$mark1[final$mark1=="SCR"] <- "SCR"
  final$mark1[final$mark1=="TRS"] <- "TRS"

  final$mark2[final$mark2=="AR"] <- "AR"
  final$mark2[final$mark2=="ARC"] <- "AR"
  final$mark2[final$mark2=="CTCF"] <- "CTCF"
  final$mark2[final$mark2=="CTCFC"] <- "CTCF"
  final$mark2[final$mark2=="EAR"] <- "EAR"
  final$mark2[final$mark2=="EARC"] <- "EAR"
  final$mark2[final$mark2=="EPR"] <- "EPR"
  final$mark2[final$mark2=="EPRC"] <- "EPR"
  final$mark2[final$mark2=="ER"] <- "ER"
  final$mark2[final$mark2=="ERC"] <- "ER"
  final$mark2[final$mark2=="HET"] <- "HET"
  final$mark2[final$mark2=="PAR"] <- "PAR"
  final$mark2[final$mark2=="PARC"] <- "PAR"
  final$mark2[final$mark2=="PPR"] <- "PPR"
  final$mark2[final$mark2=="PPRC"] <- "PPR"
  final$mark2[final$mark2=="PR"] <- "PR"
  final$mark2[final$mark2=="PRC"] <- "PR"
  final$mark2[final$mark2=="RPS"] <- "RPS"
  final$mark2[final$mark2=="SCR"] <- "SCR"
  final$mark2[final$mark2=="TRS"] <- "TRS"

  final <- GRanges(seqname = final$region1chrom,
                   ranges = IRanges(start = final$region1start,
                                    end = final$region1end),
                   mark1 = final$mark1,
                   gene1 = final$gene1,
                   FPKM1 = final$FPKM1,
                   region2chrom = final$region2chrom,
                   region2start = final$region2start,
                   region2end = final$region2end,
                   mark2 = final$mark2,
                   gene2 = final$gene2,
                   FPKM2 = final$FPKM2,
                   score = final$score)
  final
}
