% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_generic.R
\name{reduce_ggplot_object_size}
\alias{reduce_ggplot_object_size}
\title{Reduce the memory footprint of ggplot objects}
\usage{
reduce_ggplot_object_size(x)
}
\arguments{
\item{x}{a ggplot object, or a (nested) list of ggplot objects}
}
\description{
it seems like ggplot objects retain copies of source data, which can result in out of control RAM and RData file sizes.
For example, Cscore histograms from 200 Spectronaut files (i.e. 200 basic histograms of ~50k peptides) will results in 100+ GB RData files
}
\details{
Reference; https://stackoverflow.com/questions/75594274/why-are-my-ggplot-objects-getting-so-large-in-size
Reference; https://github.com/tidyverse/ggplot2/issues/4056#issue-636077598
}
