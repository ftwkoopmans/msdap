
#' returns all normalization algorithms integrated with MS-DAP
#'
#' @description
#'
#' median: scale each sample such that median abundance values are the same for all samples in the dataset.
#'
#' loess: Loess normalization as implemented in the limma R package (PMID:25605792) <https://bioconductor.org/packages/release/bioc/html/limma.html>. code: `limma::normalizeCyclicLoess(log2_data, iterations = 10, method = "fast")`. Normalize the columns of a matrix, cyclicly applying loess normalization to normalize each columns against the average over all columns.
#'
#' vsn: Variance Stabilizing Normalization (VSN) as implemented in the vsn R package (PMID:12169536) <https://bioconductor.org/packages/release/bioc/html/vsn.html>. code: `vsn::justvsn()`. From bioconductor: "The model incorporates data calibration step (a.k.a. normalization), a model for the dependence of the variance on the mean intensity and a variance stabilizing data transformation. Differences between transformed intensities are analogous to 'normalized log-ratios'".
#'
#' rlr: Robust Linear Regression normalization, as implemented in the MSqRob package (PMID:26566788) <https://github.com/statOmics/msqrob>. For each sample s, perform a robust linear regression of all values (peptide intensities) against overall median values (e.g. median value of each peptide over all samples) to obtain the normalization factor for sample s.
#'
#' msempire: log-foldchange mode normalization, as implemented in the msEmpiRe package (PMID:31235637) <https://github.com/zimmerlab/MS-EmpiRe>. Instead of computing all pairwise sample scaling factors (i.e. foldchange distributions between all pairs of samples within a sample group), MS-EmpiRe uses single linkage clustering to normalize to subsets of 'most similar' samples and iteratively expands until all within-group samples are covered.
#'
#' vwmb: Variation Within, Mode Between (VWMB) normalization. In brief, this minimizes the median peptide variation within each sample group, then scales between all pairwise sample groups such that the log-foldchange mode is zero.
#' The normalization algorithm consists of two consecutive steps:
#' 1) samples are scaled within each group such that the median of variation estimates for all rows is minimized
#' 2) summarize all samples per group by respective row mean values (from `row*sample` to a `row*group` matrix). Then rescale at the sample-group-level to minimize the mode log-foldchange between all groups
#' See further MS-DAP function \code{\link{normalize_vwmb}}.
#'
#' mwmb: Mode Within, Mode Between (MWMB) normalization. A variant of VWMB. Normalize (/scale) samples within each sample group such that their pairwise log-foldchange modes are zero, then scales between groups such that the log-foldchange mode is zero (i.e. the between-group part is the same as VWMB). If the dataset has (unknown) covariates and a sufficient number of replicates, this might be beneficial because covariate-specific effects are not averaged out as they might be with `VWMB`. See further MS-DAP function \code{\link{normalize_vwmb}}.
#'
#' modebetween: only the "Mode Between" part of VWMB described earlier, does not affect scaling between (replicate) samples within the same sample group. Note that this is mode-between normalization at peptide level, in most cases you'll want to use "modebetween_protein" instead.
#'
#' modebetween_protein  (also referred to as "MBprot", e.g. in the MS-DAP manuscript and some documentation): only the "Mode Between" part of VWMB described earlier, but the scaling factors are computed at protein-level !!  When this normalization function is used, the \code{\link{normalize_modebetween_protein}} function will first rollup the peptide data matrix to protein-level, then compute between-sample-group scaling factors and finally apply those to the input peptide-level data matrix to compute the normalized peptide data.
#'
#' vw: only perform peptide-level variation-within normalization. This rescales samples such that the peptide standard deviations (intensity values across samples in the same sample group) are minimized.  i.e. this is only the first part of the VWMB algorithm. Note that in most cases / experimental designs, you want to end your normalizations with modebetween_protein  (only normalization variation within group is often not the most suitable normalization solution)
#'
#' mw: analogous to the "vw" option, but here performing peptide-level normalization that uses mode normalization within groups (first part of MWMB). No between group rescaling. Same recommendations/notes as with "vm" apply
#'
#' @export
normalization_algorithms = function() {
  c("median", "loess", "vsn", "rlr", "msempire", "vwmb", "mwmb", "modebetween", "modebetween_protein", "vw", "mw")
}



#' Normalize a numeric matrix; in MS-DAP this is typically a peptide-level log2 intensity matrix
#'
#' wrapper function around various normalization algorithms. In the results, log2 intensity values below 1 are thresholded at 1
#'
#' @param x_as_log2 the numeric matrix, must be log2 transformed
#' @param algorithm the normalization algorithm to apply. Either one of the built-in options, see \code{\link{normalization_algorithms}} function focumentation, or the name of your custom normalization function (see GitHub documentation on custom normalization functions for more details)
#' @param group_by_cols a character array of the same length as the number of columns in the provided data matrix, describing the grouping of each column/sample
#' @param group_by_rows optionally, an array of grouping IDs for the rows. e.g. if input matrix is a peptide table these can be strings with protein-group identifiers. Required parameter for modebetween_protein normalization
#' @param rollup_algorithm optionally, the algorithm for combining peptides to proteins as used in normalization algorithms that require a priori rollup from peptides to a protein-level abundance matrix (e.g. modebetween_protein). Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#'
#' @importFrom limma normalizeCyclicLoess
#' @importFrom vsn justvsn
#' @seealso `normalization_algorithms()` for available normalization algorithms and documentation.
#' @export
normalize_matrix = function(x_as_log2, algorithm, group_by_cols = NA, group_by_rows = NA, rollup_algorithm = "maxlfq") {
  #input validation; not a valid input object, silently return
  if (!is.matrix(x_as_log2) || nrow(x_as_log2) == 0 || ncol(x_as_log2) < 2 || length(algorithm) != 1 || is.na(algorithm) || algorithm == "") {
    return(x_as_log2)
  }

  if (any(rowSums(!is.na(x_as_log2)) == 0)) {
    append_log("Remove all rows that have only NA values prior to normalization", type = "error")
  }

  x_as_log2[!is.finite(x_as_log2)] = NA
  valid_mask = length(group_by_cols) == ncol(x_as_log2) && !anyNA(group_by_cols) && (all(is.character(group_by_cols) & group_by_cols != "") || all(is.integer(group_by_cols)))
  valid_group_by_rows = length(group_by_rows) == nrow(x_as_log2) && !anyNA(group_by_rows) && (all(is.character(group_by_rows) & group_by_rows != "") || all(is.integer(group_by_rows)))

  if (algorithm == "median") {
    return(threshold_numerics(normalize_median(x_as_log2), 1))
  }

  if (algorithm == "loess") {
    # limma package
    return(threshold_numerics(limma::normalizeCyclicLoess(x_as_log2, iterations = 10, method = "fast"), 1))
  }

  if (algorithm == "rlr") {
    # implementation from MSqRob, application of MASS package
    return(threshold_numerics(normalize_rlr(x_as_log2), 1))
  }

  if (algorithm == "vsn") {
    # vsn package
    return(threshold_numerics(suppressMessages(vsn::justvsn(2^x_as_log2)), 1)) # justvsn wants non-log input
  }

  if (algorithm == "vwmb") {
    if (!valid_mask) {
      append_log("vwmb normalization requires a definition of sample groups for within- and between-group normalization", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups = group_by_cols, metric_within = "var", metric_between = "mode"), 1))
  }

  if (algorithm == "mwmb") {
    if (!valid_mask) {
      append_log("mwmb normalization requires a definition of sample groups for within- and between-group normalization", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups = group_by_cols, metric_within = "mode", metric_between = "mode"), 1))
  }

  # only normalize between replicates by minimizing variation
  if (algorithm == "vw") {
    if (!valid_mask) {
      append_log("vw (variation-within) normalization requires a definition of sample groups", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups = group_by_cols, metric_within = "var", metric_between = ""), 1))
  }

  # only normalize between replicates by minimizing mode-foldchange
  if (algorithm == "mw") {
    if (!valid_mask) {
      append_log("mw (mode-within) normalization requires a definition of sample groups", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups = group_by_cols, metric_within = "mode", metric_between = ""), 1))
  }

  # peptide-level variant of modebetween
  if (algorithm == "modebetween") {
    if (!valid_mask) {
      append_log("'mode between' normalization requires a definition of sample groups for between-group normalization", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups = group_by_cols, metric_within = "", metric_between = "mode"), 1)) # disable within-group normalization
  }

  if (algorithm == "modebetween_protein") {
    if (!valid_mask) {
      append_log("'modebetween_protein' normalization requires a definition of sample groups for between-group normalization", type = "error")
    }
    if (!valid_group_by_rows) {
      append_log("'modebetween_protein' normalization requires an array of protein identifiers (not NA or empty string), equal to the number of rows of the peptide matrix, that will be used for peptide-to-protein rollup", type = "error")
    }

    return(threshold_numerics( normalize_modebetween_protein(x_as_log2, proteins = group_by_rows, groups = group_by_cols, rollup_algorithm = rollup_algorithm), 1))
  }

  if (algorithm == "msempire") {
    if (!valid_mask) {
      append_log("MS-EmpiRe normalization requires a definition of sample groups for within-group normalization", type = "error")
    }
    if (length(unique(group_by_cols)) != 2) {
      append_log("MS-EmpiRe normalization only supports normalization between 2 groups", type = "error")
    }
    capture.output(result <- log2(suppressMessages(normalize_msempire(2^x_as_log2, group_by_cols))))
    return(threshold_numerics(result, 1))
  }

  # guess if the user provided some custom function for normalization
  f = tryCatch(match.fun(algorithm, descend = FALSE), error = function(...) NULL)
  if(!is.function(f)) {
    if(grepl("::", algorithm, fixed = T)) {
      f = tryCatch(utils::getFromNamespace(gsub(".*::", "", algorithm), gsub("::.*", "", algorithm), envir = .GlobalEnv), error = function(...) NULL)
    } else {
      f = tryCatch(utils::getFromNamespace(algorithm, "msdap", envir = .GlobalEnv), error = function(...) NULL)
    }
  }
  if(is.function(f)) {
    result = f(x_as_log2 = x_as_log2, group_by_cols = group_by_cols, group_by_rows = group_by_rows, rollup_algorithm = rollup_algorithm)
    # validation checks on expected output, to facilitate debugging/feedback for custom implementations
    if(!is.matrix(result) || nrow(result) != nrow(x_as_log2) || ncol(result) != ncol(x_as_log2) || any(colnames(result) != colnames(x_as_log2)) || mode(result) != "numeric") {
      append_log(sprintf("provided custom function for normalization '%s' must return a numeric matrix (either NA or numeric values, no infinites) with the same dimensions and column names as the input matrix", algorithm), type = "error")
    }
    return(threshold_numerics(result, 1))
  }

  # fall-through for unknown params
  append_log(paste("unsupported normalization parameter:", algorithm), type = "error")
}



#' scale each sample such that median abundance values are the same for all samples in the dataset
#'
#' @param x_as_log2 log2 transformed abundance values
#' @param ... remaining parameters are ignored
#' @importFrom matrixStats colMedians
normalize_median = function(x_as_log2, ...) {
  s = matrixStats::colMedians(x_as_log2, na.rm=T)
  s = s - mean(s)
  # slightly faster than looping columns  or  t(t(x_as_log2) - scale_per_sample)
  return(x_as_log2 - rep(s, rep.int(nrow(x_as_log2), ncol(x_as_log2))) )
}



#' wrapper function for msEmpiRe::normalize()
#'
#' @param x numerical matrix, NOT log transformed
#' @param group_by_cols a character array of the same length as the number of columns in the provided data matrix, describing the grouping of each column/sample
#'
#' @importFrom Biobase ExpressionSet annotatedDataFrameFrom fData pData exprs
#' @importFrom msEmpiRe normalize
normalize_msempire = function(x, group_by_cols) {
  # create mock Biobase::ExpressionSet
  eset_input = Biobase::ExpressionSet(
    assayData = x,
    featureData = Biobase::annotatedDataFrameFrom(x, byrow = T),
    protocolData = Biobase::annotatedDataFrameFrom(x, byrow = F)
  )
  rn = paste0("s", 1:nrow(x))
  Biobase::fData(eset_input) = data.frame(sequence = rn, row.names = rn, stringsAsFactors = T)
  Biobase::pData(eset_input) = data.frame(samples = colnames(x), condition = group_by_cols, row.names = colnames(x), stringsAsFactors = F)

  # actual normalization
  eset_norm = msEmpiRe::normalize(eset_input)
  x_norm = Biobase::exprs(eset_norm)
  x_norm[!is.finite(x_norm) | x_norm == 0] = NA
  return(x_norm)
}



#' Robust Linear Regression (RLR) normalization as sourced from MSqRob package
#'
#' copy/paste of code @ https://github.com/statOmics/MSqRob/blob/MSqRob0.7.6/R/preprocess_MaxQuant.R
#'
#' For each sample s, perform a robust linear regression of all values (peptide intensities) against overall median values (e.g. median value of each peptide over all samples) to obtain the normalization factor for sample s.
#'
#' @param exprs numerical matrix
#' @param weights optional vector of weights (per row) to be supplied to MASS::rlm()
#'
#' @importFrom MASS rlm
normalize_rlr_MSqRob_implementation = function(exprs, weights = NULL) {
  mediandata = apply(exprs, 1, stats::median, na.rm = TRUE)
  flag1 = 1
  for (j in 1:ncol(exprs)) {
    LRfit = MASS::rlm(as.matrix(exprs[, j]) ~ mediandata, weights = weights, na.action = na.exclude)
    Coeffs = LRfit$coefficients
    a = Coeffs[2]
    b = Coeffs[1]
    if (flag1 == 1) {
      globalfittedRLR = (exprs[, j] - b) / a
      flag1 = 2
    }
    else {
      globalfittedRLR = cbind(globalfittedRLR, (exprs[, j] - b) / a)
    }
  }
  return(globalfittedRLR)
}



#' Wrapper function for normalize_rlr_MSqRob_implementation
#'
#' @param x log2 abundance matrix
normalize_rlr = function(x) {
  y = normalize_rlr_MSqRob_implementation(x)
  dimnames(y) = dimnames(x)
  return(y)
}
