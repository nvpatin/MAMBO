% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlierLoadings.R
\name{outlierLoadings}
\alias{outlierLoadings}
\title{Outlier principal component loadings}
\usage{
outlierLoadings(
  results,
  locus,
  pc,
  type = c("z", "iqr"),
  thresh = 3,
  min.pct.reps = 0.95
)
}
\arguments{
\item{results}{output of a \link{mambo} run.}

\item{locus}{label name of response or predictor locus.}

\item{pc}{principal component to summarize (1:n).#'}

\item{type}{choose outlier based on inter-quartile interval (\code{iqr}) or
z-score (\code{z}).}

\item{thresh}{if \code{type = 'z'}, the minimum z-score used to identify
outliers.}

\item{min.pct.reps}{minimum percent of replicates that an ASV should be 
identified as an outlier in to be saved.}
}
\description{
Return loading summary of outlier ASVs.
}
\author{
Eric Archer \email{eric.archer@noaa.gov}
}
