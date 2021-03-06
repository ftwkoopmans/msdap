% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dea.R
\name{dea}
\alias{dea}
\title{Differential expression analysis}
\usage{
dea(
  dataset,
  qval_signif = 0.01,
  fc_signif = 0,
  algo_de = "deqms",
  algo_rollup = "maxlfq",
  output_dir_for_eset = ""
)
}
\arguments{
\item{dataset}{your dataset}

\item{qval_signif}{threshold for significance of adjusted p-values}

\item{fc_signif}{threshold for significance of log2 foldchanges. Set to zero to disregard or a positive value to apply a cutoff to absolute log2 foldchanges. MS-DAP can also perform a bootstrap analyses to infer a reasonable threshold by setting this parameter to NA}

\item{algo_de}{algorithms for differential expression analysis. options: ebayes, deqms, msqrobsum, msempire, msqrob (to run multiple, provide an array)}

\item{algo_rollup}{strategy for combining peptides to proteins as used in DEA algorithms that first combine peptides to proteins and then apply statistics, like ebayes and deqms. options: maxlfq, sum. The former applies the MaxLFQ algorithm, the latter employs the 'classic' strategy of summing all peptides per protein. See further rollup_pep2prot()}

\item{output_dir_for_eset}{optionally, provide an output directory where the expressionset objects should be stored. Only useful if you're doing downstream data analysis that requires this data}
}
\description{
Differential expression analysis
}
