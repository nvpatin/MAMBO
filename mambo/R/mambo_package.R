#' \code{mambo} package
#'
#' Metabarcoding Analysis using Modeled Bayesian Occurrences
#'
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
#' @name counts.12s
#' @title 12S ASV Read Counts
#' @description << add a description >>
#' @usage data(counts.12s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

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
#' @name counts.COI
#' @title COI ASV Read Counts
#' @description << add a description >>
#' @usage data(counts.COI)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.12s
#' @title 12S ASV Taxonomy
#' @description << add a description >>
#' @usage data(taxa.12s)
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

#' @docType data
#' @name taxa.COI
#' @title COI ASV Taxonomy
#' @description << add a description >>
#' @usage data(taxa.COI)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL
