#' \code{mambo} package
#'
#' Metabarcoding Analysis using Modeled Bayesian Occurrences
#'
#' @aliases mambo-package
#' @docType package
#' @name mambo_package
#'
#' @keywords package
#' 
#' @export
#' 
mamboTutorial <- function() {
  utils::browseURL(system.file("mambo_Tutorial.html", package = "mambo"))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Welcome to mambo v", utils::packageVersion("mambo"), "\n",
    "See mamboTutorial() for a guide to the package."
  )
}

#' @docType data
#' @name counts.16s
#' @title 16S ASV Read Counts
#' @description << add a description >>
#' @usage data(counts.16s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

#' @docType data
#' @name counts.18s
#' @title 18S ASV Read Counts
#' @description << add a description >>
#' @usage data(counts.18s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.16s
#' @title 16S ASV Taxonomy
#' @description << add a description >>
#' @usage data(taxa.16s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.18s
#' @title 18S ASV Taxonomy
#' @description << add a description >>
#' @usage data(taxa.18s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL
