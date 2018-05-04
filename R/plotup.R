#' Generates UpSetR-style plot from overlap data
#'
#' @param datafile The object output from the overlap function
#' @param ... Additional arguments passed to upset function from UpSetR package
#' @importFrom UpSetR fromList upset
#' @import utils
#' @export
#' @return Returns an UpSetR plot from overlap data
#' @examples
#' hic_chr20 <- system.file("extdata", "hic_chr20.txt", package = "HiCAGE")
#' segment_chr20 <- system.file("extdata", "segment_chr20.bed",
#'     package = "HiCAGE")
#' rna_chr20 <- system.file("extdata", "rna_chr20.tsv", package = "HiCAGE")
#' overlapoutput <- overlap(hicfile = hic_chr20,
#'                          segmentfile = segment_chr20,
#'                          rnafile = rna_chr20)
#' plotup(datafile = overlapoutput)

plotup <- function(datafile, ...) {

  m3 <- unique(append(datafile$mark1, datafile$mark2))
  marklist <- list()

  for(i in m3) {
    m1 <- which(datafile$mark1 == i)
    m2 <- which(datafile$mark2 == i)
    m3 <- c(unique(m1), unique(m2))
    m3 <- list(m3)
    names(m3) <- i
    marklist <- c(marklist, m3)
  }

  markout <- fromList(marklist)
  upset(markout, nsets = length(marklist), ...)
}
