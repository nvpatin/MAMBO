#' @title Filter count data
#' @description Remove ASVs or samples with low frequencies or read counts.
#'
#' @param x \code{data.frame} or matrix of read counts where samples are columns
#'   and rows are ASVs.
#' @param asvs.required minimum number of ASVs per sample for retention of sample
#' @param samples.required minimum number of samples an ASV occurs in for retention of ASV.
#' @param min.sample.read.count minimum total read count of a sample for retention.
#' @param min.asv.read.pct minimum percent of all reads for an ASV (across all samples) for
#'   retention of ASV.
#'
#' @return matrix of count data with filters applied.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
filterCountData <- function(x, asvs.required = 2, samples.required = 2,
                            min.sample.read.count = 1000,
                            min.asv.read.pct = 1.5e-6) {
  
  message('Input: ', ncol(x), ' samples and ', nrow(x), ' ASVs.')
  new.x <- as.matrix(x)
  
  # iteratively check that there are enough ASVs per sample and enough samples per ASV
  has.enough <- FALSE
  while(!has.enough) {
    has.enough.asvs <- apply(x, 2, function(sample.i) sum(sample.i > 0) >= asvs.required)
    new.x <- new.x[, has.enough.asvs]
    has.enough.samples <- apply(new.x, 1, function(asv.i) sum(asv.i > 0) >= samples.required)
    new.x <- new.x[has.enough.samples, ]
    has.enough <- all(has.enough.asvs) & all(has.enough.samples)
  }
  message(ncol(x) - ncol(new.x), ' samples removed due to insufficent ASVs per sample.')
  message(nrow(x) - nrow(new.x), ' ASVs removed due to insufficient samples per ASV.')
  
  # check that samples don't have too few reads overall
  sample.has.enough.reads <- colSums(x) >= min.sample.read.count
  new.x <- new.x[, sample.has.enough.reads]
  message(sum(!sample.has.enough.reads), ' samples removed due to low read counts.')
  
  # check that ASVs have enough read counts across samples
  asv.has.enough.reads <- rowSums(new.x) > (sum(new.x) * min.asv.read.pct)
  new.x <- new.x[asv.has.enough.reads, ]
  message(sum(!asv.has.enough.reads), ' ASVs removed due to low read counts across samples.')
  
  message('Retained: ', ncol(new.x), ' samples and ', nrow(new.x), ' ASVs.')
  new.x
}
