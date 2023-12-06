#' @title Load example data
#' @description Load named example data files included in package
#'
#' @param file name of data file.
#'
#' @return an ASV table as a data frame.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x <- exampleData("Merged2018_16S_otu_filtered.csv")
#' 
#' str(x)
#'
#' @export
#'
exampleData <- function(file) {
  f <- system.file(file.path("extdata", file), package = "mambo")
  if(f == "") stop("Filename '", file, "' not found.")
  read.csv(f, row.names = 1)
}
