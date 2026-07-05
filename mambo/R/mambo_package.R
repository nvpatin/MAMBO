"_PACKAGE"
#' @name mambo_package
#' @title Metabarcoding Analysis using Modeled Bayesian Occurrences
#' 
#' @keywords internal
#' @importFrom rlang .data
#' 
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Welcome to mambo v", utils::packageVersion("mambo"), "\n",
    "See mamboTutorial() for a guide to the package."
  )
}

#' @docType data
#' @name counts.12s
#' @title 12S ASV Read Counts
#' @description Blackman et al. 2022 12S ASV table
#' @usage data(counts.12s)
#' @format data.frame
#' @references https://doi.org/10.1038/s42003-022-03216-z
#' @keywords datasets
NULL

#' @docType data
#' @name counts.16s
#' @title 16S ASV Read Counts
#' @description CN18S, CN18F, and Lasker 2018 16S ASV table
#' @usage data(counts.16s)
#' @format data.frame
#' @keywords datasets
NULL

#' @docType data
#' @name counts.18s
#' @title 18S ASV Read Counts
#' @description CN18S, CN18F, and Lasker 2018 18S ASV table
#' @usage data(counts.18s)
#' @format data.frame
#' @keywords datasets
NULL

#' @docType data
#' @name counts.COI
#' @title COI ASV Read Counts
#' @description Blackman et al. 2022 COI ASV table
#' @usage data(counts.COI)
#' @format data.frame
#' @references https://doi.org/10.1038/s42003-022-03216-z
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.12s
#' @title 12S ASV Taxonomy
#' @description Blackman et al. 2022 12S taxonomy
#' @usage data(taxa.12s)
#' @format data.frame
#' @references https://doi.org/10.1038/s42003-022-03216-z
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.16s
#' @title 16S ASV Taxonomy
#' @description CN18S, CN18F, and Lasker 2018 16S taxonomy
#' @usage data(taxa.16s)
#' @format data.frame
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.18s
#' @title 18S ASV Taxonomy
#' @description CN18S, CN18F, and Lasker 2018 18S taxonomy
#' @usage data(taxa.18s)
#' @format data.frame
#' @keywords datasets
NULL

#' @docType data
#' @name taxa.COI
#' @title COI ASV Taxonomy
#' @description Blackman et al. 2022 COI taxonomy
#' @usage data(taxa.COI)
#' @format data.frame
#' @references https://doi.org/10.1038/s42003-022-03216-z
#' @keywords datasets
NULL

#' @docType data
#' @name sample.meta
#' @title Metadata for CCE samples
#' @description CN18S, CN18F, and Lasker 2018 environmental metadata
#' @usage data(sample.meta)
#' @format data.frame
#' @keywords datasets
NULL