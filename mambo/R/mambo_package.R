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
#' @name fl16s
#' @title Flyer 16S Read Counts
#' @description << add a description >>
#' @usage data(fl16s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL

#' @docType data
#' @name fl18s
#' @title Flyer 18S Read Counts
#' @description << add a description >>
#' @usage data(fl18s)
#' @format data.frame
#' @references << add a reference >>
#' @keywords datasets
NULL