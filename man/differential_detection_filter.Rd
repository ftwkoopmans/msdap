% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_differential_detect.R
\name{differential_detection_filter}
\alias{differential_detection_filter}
\title{Protein differential detection filtering using a hardcoded number of peptide observations per condition}
\usage{
differential_detection_filter(
  dataset,
  k_diff = NA,
  frac_diff = NA,
  npep_pass = 2L,
  nobs_ratio = 3,
  int_ratio = 0,
  normalize_intensities = TRUE
)
}
\arguments{
\item{dataset}{your dataset. At least 1 contrast should have been specified prior}

\item{k_diff}{peptide-level criterium: peptide must be detected in at least k more samples in condition 1 versus condition 2, or vice versa}

\item{frac_diff}{peptide-level criterium: peptide must be detected in at least x\% more samples in condition 1 versus condition 2, or vice versa}

\item{npep_pass}{minimum number of peptides that must pass the peptide-level differential detection criteria. Default: 2}

\item{nobs_ratio}{minimum enrichment ratio between experimental conditions for the total peptide detection count per protein (see output table specification, \code{nobs1} and \code{log2fc_nobs2/nobs1} to which this filter is applied). Note that this value is NOT on log2 scale, i.e. set 2 for 2-fold enrichment. Default: 4}

\item{int_ratio}{analogous to \code{nobs_ratio}, but for the enrichment in sum peptide intensity values. Default: 0 (disabled)}

\item{normalize_intensities}{normalize the protein intensity matrix prior to computing sum intensities and intensity ratios. Default: TRUE}
}
\value{
A tibble where each row describes 1 proteingroup ("protein_id") in 1 contrast, with the following columns:
\itemize{
\item \code{npep_total} = total number of peptides detected across any of the samples in the current contrast. Useful for post-hoc filtering, e.g. when you do not set stringent criteria for \code{npep_pass}
\item \code{npep_pass1} = number of peptides that pass the specified filtering rules in condition 1 of the current contrast
\item \code{npep_pass2} = analogous to \code{npep_pass1}, but for condition 2
\item \code{nobs1} = sum of peptide detections across all samples in condition 1 (i.e. each detected peptide is counted once per sample, this is the total sum across peptides*samples for respective protein_id and condition). Useful for post-hoc filtering, e.g. when you do not set stringent criteria for \code{npep_pass}
\item \code{nobs2} = analogous to \code{nobs1}, but for condition 2
\item \code{fracobs1} = percentage of all possible detections made within condition 1 = number of observed datapoints / (#peptide * #sample)
\item \code{fracobs2} = analogous to \code{fracobs1}, but for condition 2
\item \code{log2fc_nobs2/nobs1} = log2 foldchange of observation counts. Positive values are enriched in condition 2. Proteins exclusive to condition 2 have value \code{Inf} and exclusive to condition 1 have value \code{-Inf}
\item \code{int1} = sum peptide intensity across all samples in condition 1
\item \code{int2} = sum peptide intensity across all samples in condition 2
\item \code{log2fc_int2/int1} = log2 foldchange of protein intensities, analogous to \code{log2fc_nobs2/nobs1}
\item \code{pass} = protein matches input filtering
}
}
\description{
\subsection{Pseudocode}{

Number of peptides differentially detected in condition 1 as compared to condition 2:

\code{npep_pass1 = sum(n1 - n2 >= k1_diff)}

Where n1 and n2 are vectors with the number of detects per peptide in conditions 1 and 2, respectively,
and k1_diff is the sample count threshold for differential detection.

Now we can test at protein-level, directionally!
We're looking for proteins that are near-exclusively detected in condition 1 and their nobs ratio is in the same direction.

\code{npep_pass1 >= npep_pass & nobs1 > nobs_ratio * nobs2}

Where npep_pass and nobs_ratio are user-provided thresholds, nobs1 is the sum of n1 (see further return value specification for this function).
Above rule tests enrichment in condition 1, afterwards we also test analogously with conditions 1 and 2 swapped.
}

\subsection{real-world example data:}{

2 conditions, 6 replicates each:

\tabular{llrr}{
\strong{protein_id} \tab \strong{peptide_id} \tab \strong{n1} \tab \strong{n2} \cr
P09041 \tab ALMDEVVK_2 \tab 2 \tab 6 \cr
P09041 \tab LTLDKVDLK_2 \tab 1 \tab 5
}

n2-n1 diff is 4 for both peptides, nobs ratio was 3.7-fold enrichment in condition 2.
}
}
\examples{
\dontrun{
## example 1:
# default / stringent: find proteins that have at least 2 peptides
# detected in 66\% more samples in either condition (with a minimum of 3 samples)
# AND the overall detection rate is larger than 3-fold
x = differential_detection_filter(
  dataset, k_diff = 3, frac_diff = 0.66, npep_pass = 2, nobs_ratio = 3
)

# code snippet to create a pretty-print table with resulting proteins
y = x \%>\%
  # only retain proteins that match input criteria
  filter(pass) \%>\%
  # add protein metadata and rearrange columns
  left_join(dataset$proteins \%>\%
    select(protein_id, fasta_headers, gene_symbols_or_id),
    by = "protein_id") \%>\%
  select(contrast, protein_id, fasta_headers, gene_symbols_or_id,
         tidyselect::everything()) \%>\%
  # for prettyprint, trim the contrast names
  mutate(contrast = gsub(" *#.*", "", sub("^contrast: ", "", contrast))) \%>\%
  # sort data by column, then by ratio
  arrange(contrast, `log2fc_nobs2/nobs1`)
print(y)

## example 2:
# a more lenient filter: apply criteria to only 1 (or 0) peptides but
# rely mostly on the nobs_ratio criterium and additionally
# add post-hoc filtering on the total number of peptides per protein
# and require either protein to have an overall 50\% detection rate
# (across peptides and samples)
x = differential_detection_filter(
  dataset, nobs_ratio = 3, npep_pass = 1, k_diff = 3, frac_diff = 0.66
) \%>\%
  mutate(pass = pass & npep_total > 1 & pmax(fracobs1, fracobs1) >= 0.5)
}
}
