% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_vwmb.R
\name{pairwise_modes}
\alias{pairwise_modes}
\title{compute the mode of all sample/column pairs
column pairs that have less than 10 values in common are skipped (result value left at default 0)}
\usage{
pairwise_modes(x)
}
\arguments{
\item{x}{log-transformed data matrix, where columns are samples and rows are features}
}
\value{
a matrix of ncol*ncol pairwise foldchange-modes
}
\description{
compute the mode of all sample/column pairs
column pairs that have less than 10 values in common are skipped (result value left at default 0)
}
