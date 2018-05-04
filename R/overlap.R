#' Combine chromosome conformation capture data with
#' segmentation and RNA data
#'
#' @param hicfile The chromosome conformation capture datafile
#' @param segmentfile The segmentation datafile from StateHub for example
#' @param rnafile The RNA gene expression datafile containing Ensembl IDs and
#' FPKM or TPM normalized expression values
#' @param bio_mart BiomaRt used for pulling gene name and location data
#' @param martset The specific mart or dataset used for genome information
#' @param webhost Website host to connect to for downloading biomaRt data
#' @param hic.columns Columns from the chromsome conformation capture data file
#' that contain, in order, 'left chromosome', 'left start', 'left end', 'right
#' chromosome', 'right start', 'right end', 'interaction score'
#' @param segment.columns Columns from the segmentation data file that contain,
#' in order, 'Chromosome', 'Start', 'End', 'Mark', 'Score'. Select only
#' 'Chromosome', 'Start', 'End', and 'Mark' if manual prioritization is
#' desired
#' @param rna.columns Columns from the RNA-seq data file that contain, in order,
#' 'GeneID', 'FPKM'
#' @param manual.priority A concatenated list of state calls contained in the
#' segmentation file that are to be used to manually prioritize state calls that
#' fall within the same genomic region of the chromosome conformation capture
#' file. This argument can be used if the segmentation file does not contain
#' scores for state calls or if scores want to be overridden. The concatenated
#' list should be in the form c("PAR", "EAR", "AR", "PPR", "EPR", "TRS", "HET")
#' ordered from highest priority to lowest priority. In this case, only select
#' columns containing'Chromosome', 'Start', 'End', and 'Mark' in the
#' segmentation file using `segment.columns`.
#' @param prune.priority A concatenated list of state calls contained in the
#' segmentation file that are to be used to prioritize state calls that both
#' fall within the same genomic region of the chromosome conformation capture
#' file and have identical scores. In this instance, segmentation scores will
#' first be used to prioritize states, then prune.priority will be used to break
#' ties. This argument is used so that only one interaction is defined for each
#' line of the chromosome conformation capture file. Otherwise, multiple
#' interactions can potentially be defined for a single chromsome interaction if
#' two or more state calls in that region have identical segmentation scores.
#' The concatenated list should be in the form c("PAR", "EAR", "AR", "PPR",
#' "EPR", "TRS", "HET") ordered from highest priority to lowest priority. If
#' manual.priority is set, this argument will have no effect
#' @importFrom readr read_tsv read_csv
#' @importFrom tidyr separate
#' @importFrom dplyr mutate left_join full_join group_by bind_rows bind_cols
#' anti_join filter
#' @importFrom GenomicRanges GRanges nearest findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom biomaRt useEnsembl getBM useMart
#' @importFrom magrittr %>%
#' @export
#' @return Returns a GRanges with each chromsome conformation capture segment
#' annotated with a prioritized segmentation mark, the nearest gene to the
#' segmentation mark and the gene's FPKM expression data
#' @examples
#' hic_chr20 <- system.file("extdata", "hic_chr20.txt", package = "HiCAGE")
#' segment_chr20 <- system.file("extdata", "segment_chr20.bed",
#'     package = "HiCAGE")
#' rna_chr20 <- system.file("extdata", "rna_chr20.tsv", package = "HiCAGE")
#' example <- overlap(hicfile = hic_chr20,
#'                    segmentfile = segment_chr20,
#'                    rnafile = rna_chr20)
overlap <- function(hicfile,
                    segmentfile,
                    rnafile,
                    bio_mart = "ensembl",
                    martset = "hsapiens_gene_ensembl",
                    webhost = "www.ensembl.org",
                    hic.columns = c(1:6,8),
                    segment.columns = c(1:5),
                    rna.columns = c(1,7),
                    manual.priority = FALSE,
                    prune.priority = FALSE) {
  #Progress Bar
  pb   <- txtProgressBar(1, 100,
                         initial = 6,
                         style=3)
  # Load biomart for assigning nearest gene and subset in GRanges
  ensembl = useMart(biomart = bio_mart,
                    dataset = martset,
                    host = webhost)
  setTxtProgressBar(pb, 31)
  gene <- getBM(attributes = c('ensembl_gene_id',
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                mart = ensembl)

  gene$chromosome_name <- gsub("^chr", "", gene$chromosome_name)
  ensGene <- GRanges(seqnames = gene$chromosome_name,
                     ranges = IRanges(start = as.integer(gene$start_position),
                                      end = as.integer(gene$end_position)))
  setTxtProgressBar(pb, 60)
  # Manual select allows user to manual choose columns from Hi-C datafile
    HiCdata <- read_tsv(file = hicfile,
                        comment = "#",
                        guess_max = 100000)
    HiCdata <- subset(HiCdata, select = hic.columns)
    colnames(HiCdata) <- c("region1chrom",
                           "region1start",
                           "region1end",
                           "region2chrom",
                           "region2start",
                           "region2end",
                           "score")
  #Remove duplicates from Hi-C data and subset into GRanges
    HiCdata <- HiCdata[!duplicated(HiCdata[c("region1chrom",
                                             "region1start",
                                             "region1end",
                                             "region2chrom",
                                             "region2start",
                                             "region2end")]),]

  HiCdata$region1chrom <- gsub("^chr", "", HiCdata$region1chrom)
  HiCdata$region2chrom <- gsub("^chr", "", HiCdata$region2chrom)

  chiaregion1 <- subset(HiCdata, select = c(1:3))
  chiaregion1 <- GRanges(seqnames = chiaregion1$region1chrom,
                         ranges =
                           IRanges(start = as.integer(chiaregion1$region1start),
                                   end = as.integer(chiaregion1$region1end)))

  chiaregion2 <- subset(HiCdata, select = c(4:6))
  chiaregion2 <- GRanges(seqnames = chiaregion2$region2chrom,
                         ranges =
                           IRanges(start = as.integer(chiaregion2$region2start),
                                   end = as.integer(chiaregion2$region2end)))

  #Parse Segmentation data and subset in GenomicRanges
  #Segmentation files from StateHub Default Model
  #http://statehub.org/modeltracks/default_model/
  segmentation <- read_tsv(file = segmentfile,
                           col_names = FALSE,
                           comment = "#",
                           skip = 1,
                           guess_max = 100000)
  segmentation <- subset(segmentation, select = segment.columns)
  if (all(manual.priority != FALSE)) {
    length(unique(segmentation$X4))
    scorelist <- list()
    for (i in 1:length(unique(segmentation$X4))) {
      score <- 1/i
      scorelist <- c(scorelist, score)
    }
    priority.list <- as.list(manual.priority)
    priority <- do.call(rbind, Map(data.frame, state = priority.list, segscore = scorelist))
    segmentation <- left_join(segmentation, priority, by = c("X4" = "state"))
  }
  colnames(segmentation) <- c("chromosome",
                              "segstart",
                              "segend",
                              "state",
                              "segscore")

  segmentation$chromosome <- gsub("^chr", "", segmentation$chromosome)
  epigenetic <- GRanges(seqnames = segmentation$chromosome,
                        ranges = IRanges(start =
                                           as.integer(segmentation$segstart),
                                         end =
                                           as.integer(segmentation$segend)),
                        mark = as.character(segmentation$state))
  if (is.null(rnafile)) {
    #Intersect ChIA-PET data with Epigenetic Segmentation Data

    overlap1 <- data.frame(findOverlaps(chiaregion1, epigenetic))

    overlap2 <- data.frame(findOverlaps(chiaregion2, epigenetic))

    region1data <- data.frame(overlap1[, "queryHits"],
                              HiCdata[overlap1$queryHits, "region1chrom"],
                              HiCdata[overlap1$queryHits, "region1start"],
                              HiCdata[overlap1$queryHits, "region1end"],
                              segmentation[overlap1$subjectHits,"state"],
                              segmentation[overlap1$subjectHits,"segscore"])
    colnames(region1data) <- c("id",
                               "region1chrom",
                               "region1start",
                               "region1end",
                               "state1",
                               "segscore1")

    region2data <- data.frame(overlap2[, "queryHits"],
                              HiCdata[overlap2$queryHits, "region2chrom"],
                              HiCdata[overlap2$queryHits, "region2start"],
                              HiCdata[overlap2$queryHits, "region2end"],
                              segmentation[overlap2$subjectHits,"state"],
                              segmentation[overlap2$subjectHits,"segscore"],
                              HiCdata[overlap2$queryHits, "score"])
    colnames(region2data) <- c("id",
                               "region2chrom",
                               "region2start",
                               "region2end",
                               "state2",
                               "segscore2",
                               "score")

    xleft <- subset(HiCdata, select = c(1:3))
    yright <- subset(HiCdata, select = c(4:7))

    xleft <- mutate(xleft, id=rownames(xleft))
    xleft$id = as.integer(xleft$id)
    yright <- mutate(yright, id=rownames(yright))
    yright$id = as.integer(yright$id)

    xleft <- anti_join(xleft, region1data, by = "id")
    yright <- anti_join(yright, region2data, by = "id")

    xleft <- mutate(xleft, state1 = "None")
    yright <- mutate(yright, state2 = "None")

    region1data$state1 <- as.character(region1data$state1)
    region1data$region1chrom <- as.character(region1data$region1chrom)
    region2data$state2 <- as.character(region2data$state2)
    region2data$region2chrom <- as.character(region2data$region2chrom)

    region1data <- bind_rows(region1data, xleft)
    region2data <- bind_rows(region2data, yright)

    region1data[is.na(region1data)] <- 0
    region2data[is.na(region2data)] <- 0

    region1data <- region1data %>%
      group_by(id) %>%
      filter(segscore1 == max(segscore1))
    region1data <- region1data[!duplicated(region1data[c("id", "state1")]),]
    region2data <- region2data %>%
      group_by(id) %>%
      filter(segscore2 == max(segscore2))
    region2data <- region2data[!duplicated(region2data[c("id", "state2")]),]


    final <- as.data.frame(full_join(region1data, region2data, by = "id"))

    final <- setNames(final, c("ID",
                               "region1chrom",
                               "region1start",
                               "region1end",
                               "mark1",
                               "segscore1",
                               "region2chrom",
                               "region2start",
                               "region2end",
                               "mark2",
                               "segscore2",
                               "HiCscore"))
  }

  else{
    #Parse RNA seq file
    RNAseq <- read_tsv(rnafile,
                       comment = "#",
                       guess_max = 100000) %>%
      subset(select = rna.columns) %>%
      separate(col = 1, into = c("gene_id", "extraint"), sep = "\\.") %>%
      subset(select = -extraint)
    colnames(RNAseq) <- c("gene_id",
                          "FPKM")
    setTxtProgressBar(pb, 75)


    #Find nearest gene to each segmentation and Hi-C location; Add FPKM
    epineargene <- data.frame(nearest(epigenetic, ensGene, ignore.strand=TRUE))
    colnames(epineargene) <- "ensembl"
    epineargene <- data.frame(gene[epineargene$ensembl, "ensembl_gene_id"])
    colnames(epineargene) <- "ensembl"
    epineargene$ensembl <- as.character(epineargene$ensembl)
    epineargene <- left_join(epineargene, RNAseq, by = c("ensembl" = "gene_id"))
    segmentation <- cbind(segmentation, epineargene)

    HiCneargene1 <- data.frame(nearest(chiaregion1, ensGene,
                                       ignore.strand=TRUE))
    colnames(HiCneargene1) <- "ensembl"
    HiCneargene1 <- data.frame(gene[HiCneargene1$ensembl, "ensembl_gene_id"])
    colnames(HiCneargene1) <- "ensembl"
    HiCneargene1$ensembl <- as.character(HiCneargene1$ensembl)
    HiCneargene1 <- left_join(HiCneargene1, RNAseq,
                              by = c("ensembl" = "gene_id"))
    setTxtProgressBar(pb, 85)

    HiCneargene2 <- data.frame(nearest(chiaregion2, ensGene,
                                       ignore.strand=TRUE))
    colnames(HiCneargene2) <- "ensembl"
    HiCneargene2 <- data.frame(gene[HiCneargene2$ensembl, "ensembl_gene_id"])
    colnames(HiCneargene2) <- "ensembl"
    HiCneargene2$ensembl <- as.character(HiCneargene2$ensembl)
    HiCneargene2 <- left_join(HiCneargene2, RNAseq,
                              by = c("ensembl" = "gene_id"))

  #Intersect ChIA-PET data with Epigenetic Segmentation Data

  overlap1 <- data.frame(findOverlaps(chiaregion1, epigenetic))

  overlap2 <- data.frame(findOverlaps(chiaregion2, epigenetic))

  region1data <- data.frame(overlap1[, "queryHits"],
                            HiCdata[overlap1$queryHits, "region1chrom"],
                            HiCdata[overlap1$queryHits, "region1start"],
                            HiCdata[overlap1$queryHits, "region1end"],
                            segmentation[overlap1$subjectHits,"state"],
                            segmentation[overlap1$subjectHits,"segscore"],
                            segmentation[overlap1$subjectHits,"ensembl"],
                            segmentation[overlap1$subjectHits, "FPKM"])
  colnames(region1data) <- c("id",
                             "region1chrom",
                             "region1start",
                             "region1end",
                             "state1",
                             "segscore1",
                             "ensembl1",
                             "FPKM1")

  region2data <- data.frame(overlap2[, "queryHits"],
                            HiCdata[overlap2$queryHits, "region2chrom"],
                            HiCdata[overlap2$queryHits, "region2start"],
                            HiCdata[overlap2$queryHits, "region2end"],
                            segmentation[overlap2$subjectHits,"state"],
                            segmentation[overlap2$subjectHits,"segscore"],
                            segmentation[overlap2$subjectHits,"ensembl"],
                            segmentation[overlap2$subjectHits, "FPKM"],
                            HiCdata[overlap2$queryHits, "score"])
  colnames(region2data) <- c("id",
                             "region2chrom",
                             "region2start",
                             "region2end",
                             "state2",
                             "segscore2",
                             "ensembl2",
                             "FPKM2",
                             "score")


  xleft <- bind_cols(HiCneargene1, HiCdata)
  xleft <- subset(xleft, select = c(1:5))
  yright <- bind_cols(HiCneargene2, HiCdata)
  yright <- subset(yright, select = c(1,2,6:9))

  xleft <- mutate(xleft, id=rownames(xleft))
  xleft$id = as.integer(xleft$id)
  yright <- mutate(yright, id=rownames(yright))
  yright$id = as.integer(yright$id)

  xleft <- anti_join(xleft, region1data, by = "id")
  yright <- anti_join(yright, region2data, by = "id")

  xleft <- mutate(xleft, state1 = "None")
  yright <- mutate(yright, state2 = "None")

  colnames(xleft)[1:2] <- c("ensembl1", "FPKM1")
  colnames(yright)[1:2] <- c("ensembl2", "FPKM2")

  region1data$ensembl1 <- as.character(region1data$ensembl1)
  region1data$state1 <- as.character(region1data$state1)
  region1data$region1chrom <- as.character(region1data$region1chrom)
  region2data$ensembl2 <- as.character(region2data$ensembl2)
  region2data$state2 <- as.character(region2data$state2)
  region2data$region2chrom <- as.character(region2data$region2chrom)

  region1data <- bind_rows(region1data, xleft)
  region2data <- bind_rows(region2data, yright)

  region1data[is.na(region1data)] <- 0
  region2data[is.na(region2data)] <- 0

  if (all(prune.priority == FALSE)) {
  region1data <- region1data %>%
    group_by(id) %>%
    filter(segscore1 == max(segscore1))
  region1data <- region1data[!duplicated(region1data[c("id", "state1")]),]
  setTxtProgressBar(pb, 90)
  region2data <- region2data %>%
    group_by(id) %>%
    filter(segscore2 == max(segscore2))
  region2data <- region2data[!duplicated(region2data[c("id", "state2")]),]
  }
  else {
    scorelist <- list()
    for (i in 1:length(unique(c(region1data$state1,region2data$state2)))) {
      score <- 1/i
      scorelist <- c(scorelist, score)
    }
    prune.list <- as.list(prune.priority)
    priority <- do.call(rbind, Map(data.frame, state = prune.list, prunescore = scorelist))
    priority$state <- as.character(priority$state)
    region1data <- left_join(region1data, priority, by = c("state1" = "state"))
    region2data <- left_join(region2data, priority, by = c("state2" = "state"))

    region1data <- region1data %>%
      group_by(id) %>%
      filter(segscore1 == max(segscore1)) %>%
      filter(prunescore == max(prunescore))
    region1data <- region1data[!duplicated(region1data[c("id", "state1")]),]
    region1data$prunescore <- NULL

    region2data <- region2data %>%
      group_by(id) %>%
      filter(segscore2 == max(segscore2)) %>%
      filter(prunescore == max(prunescore))
    region2data <- region2data[!duplicated(region2data[c("id", "state2")]),]
    region2data$prunescore <- NULL

  }

  final <- as.data.frame(full_join(region1data, region2data, by = "id"))

  final$FPKM1 <- log(final$FPKM1+1)
  final$FPKM2 <- log(final$FPKM2+1)

  final <- setNames(final, c("ID",
                             "region1chrom",
                             "region1start",
                             "region1end",
                             "mark1",
                             "segscore1",
                             "gene1",
                             "logFPKM1",
                             "region2chrom",
                             "region2start",
                             "region2end",
                             "mark2",
                             "segscore2",
                             "gene2",
                             "logFPKM2",
                             "HiCscore"))
  }

  final <- final[!duplicated(final[c("region1chrom",
                                     "region1start",
                                     "region1end",
                                     "mark1",
                                     "region2chrom",
                                     "region2start",
                                     "region2end",
                                     "mark2")]),]

  final <- final[,-1]
  setTxtProgressBar(pb, 100)
  final
}
