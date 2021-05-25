
#' the set of available normalization algorithms
#' @export
normalization_algorithms = function() {
  c("median", "loess", "vsn", "rlr", "msempire", "vwmb", "modebetween", "modebetween_protein")
}



#' Normalize a numeric matrix
#'
#' wrapper function around various normalization algorithms. Thresholds returned log2 intensity values below 1.
#'
#' @param x_as_log2 the numeric matrix, must already be log2 transformed
#' @param algorithm the normalization algorithm to apply. Either one of the built-in options, see normalization_algorithms(), or the name of your custom normalization function (see online documentation/examples)
#' @param mask_sample_groups a character array of the same length as the number of columns in the provided data matrix, describing the grouping of each column/sample
#'
#' @importFrom limma normalizeCyclicLoess
#' @importFrom vsn justvsn
#' @export
normalize_matrix = function(x_as_log2, algorithm, mask_sample_groups = NA, rows_group_by = NA) {
  #input validation; not a valid input object, silently return
  if (!is.matrix(x_as_log2) || nrow(x_as_log2) == 0 || ncol(x_as_log2) < 2 || length(algorithm) != 1 || is.na(algorithm) || algorithm == "") {
    return(x_as_log2)
  }

  if (any(rowSums(!is.na(x_as_log2)) == 0)) {
    append_log("Remove all rows that have only NA values prior to normalization", type = "error")
  }

  x_as_log2[!is.finite(x_as_log2)] = NA
  valid_mask = length(mask_sample_groups) == ncol(x_as_log2) && all(!is.na(mask_sample_groups) & mask_sample_groups != "")
  valid_rows_group_by = length(rows_group_by) == nrow(x_as_log2) && all(!is.na(rows_group_by) & rows_group_by != "")

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
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups=mask_sample_groups), 1)) # use default settings
  }

  if (algorithm == "modebetween") {
    if (!valid_mask) {
      append_log("'mode between' normalization requires a definition of sample groups for between-group normalization", type = "error")
    }
    return(threshold_numerics(normalize_vwmb(x_as_log2, groups=mask_sample_groups, metric_within = ""), 1)) # disable within-group normalization
  }

  if (algorithm == "modebetween_protein") {
    if (!valid_mask) {
      append_log("'modebetween_protein' normalization requires a definition of sample groups for between-group normalization", type = "error")
    }
    if (!valid_rows_group_by) {
      append_log("'modebetween_protein' normalization requires an array of protein identifiers (not NA or empty string), equal to the number of rows of the peptide matrix, that will be used for peptide-to-protein rollup", type = "error")
    }

    ## rollup
    # from peptide matrix to long-format tibble
    x_as_log2_tib = as_tibble(x_as_log2) %>%
      add_column(peptide_id=1:nrow(x_as_log2),
                 protein_id=rows_group_by) %>%
      tidyr::pivot_longer(cols = c(-protein_id, -peptide_id), names_to = "sample_id", values_to = "intensity") %>%
      filter(is.finite(intensity) & intensity > 0)
    # rollup by MaxLFQ
    x_as_log2__rollup = rollup_pep2prot_maxlfq(x_as_log2_tib, intensity_is_log2 = TRUE, implementation = "iq", return_as_matrix = TRUE)

    # x_as_log2__rollup = rollup_pep2prot_summation(x_as_log2_tib, intensity_is_log2 = TRUE, return_as_matrix = TRUE)
    # oneliner for summation-based rollup; x_as_log2__rollup = suppressMessages(as_tibble(2^x_as_log2, .name_repair = "universal") %>% replace(is.na(.), 0) %>% add_column(protein_id=rows_group_by) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA) %>% select(-protein_id) %>% log2 %>% as.matrix)

    # importantly, enforce consistent column order to preserve aligned with sample-to-group assignments
    stopifnot(ncol(x_as_log2) == ncol(x_as_log2__rollup) && all(colnames(x_as_log2) %in% colnames(x_as_log2__rollup)))
    x_as_log2__rollup = x_as_log2__rollup[,match(colnames(x_as_log2), colnames(x_as_log2__rollup)),drop=F] # use match() instead of direct key/string-based indexing because some samples may have names like 1,2,3,4 (eg; if key_sample is used for column names instead of sample_id, as we do in filter_dataset() )

    # vwmb + param to return scaling factors
    x_as_log2__rollup = normalize_vwmb(x_as_log2__rollup, groups=mask_sample_groups, metric_within = "", include_attributes = TRUE) # disable within-group normalization
    scale_per_sample = attr(x_as_log2__rollup, "scaling")

    # apply scaling factors, determined through normalization on protein-level, to peptide-level (simple loop is more efficient than transpose+scale+transpose)
    for(j in seq_along(scale_per_sample)) {
      x_as_log2[,j] = x_as_log2[,j] - scale_per_sample[j]
    }
    return(threshold_numerics(x_as_log2, 1))
    # debug; x_as_log2 = cbind(rnorm(100,2), rnorm(100,2.5), rnorm(100,4.5), rnorm(100,5)); mask_sample_groups = c("a","a","b","b"); rows_group_by = rep(1:(nrow(x_as_log2)/2), times=2); boxplot(x_as_log2)
  }

  if (algorithm == "msempire") {
    if (!valid_mask) {
      append_log("MS-EmpiRe normalization requires a definition of sample groups for within-group normalization", type = "error")
    }
    if (length(unique(mask_sample_groups)) != 2) {
      append_log("MS-EmpiRe normalization only supports normalization between 2 groups", type = "error")
    }
    capture.output(result <- log2(suppressMessages(normalize_msempire(2^x_as_log2, mask_sample_groups))))
    return(threshold_numerics(result, 1))
  }

  # guess if the user provided some custom function for normalization
  if (algorithm %in% ls(envir=.GlobalEnv)) {
    f = match.fun(algorithm)
    result = f(x_as_log2 = x_as_log2, mask_sample_groups = mask_sample_groups, rows_group_by = rows_group_by)
    # validation checks on expected output, to facilitate debugging/feedback for custom implementations
    if(!is.matrix(result) || nrow(result) != nrow(x_as_log2) || ncol(result) != ncol(x_as_log2) || colnames(result) != colnames(x_as_log2) || mode(result) != "numeric" || any(!is.na(result) & is.infinite(result))) {
      append_log(sprintf("provided custom function for normalization '%s' must return a numeric matrix (either NA or numeric values, no infinites) with the same dimensions and column names as the input matrix", algorithm), type = "error")
    }
    return(threshold_numerics(result, 1))
  }

  # fall-through for unknown params
  append_log(paste("unsupported normalization parameter:", algorithm), type = "error")
}



#' simply scale samples by median value
#' @param x_as_log2 log2 transformed abundance values
#' @param ... remaining parameters are ignored
#' @importFrom matrixStats colMedians
normalize_median = function(x_as_log2, ...) {
  s = matrixStats::colMedians(x_as_log2, na.rm=T)
  s = s - mean(s)
  # slightly faster than looping columns  or  t(t(x_as_log2) - scale_per_sample)
  return(x_as_log2 - rep(s, rep.int(nrow(x_as_log2), ncol(x_as_log2))) )
}



#' adapted from msEmpiRe::normalize(), supports any input matrix + definition of groups
#' input expected as finite values, not log transformed
#' @param x todo
#' @param mask_sample_groups todo
#'
#' @importFrom Biobase ExpressionSet annotatedDataFrameFrom fData pData exprs
#' @importFrom msEmpiRe normalize
normalize_msempire = function(x, mask_sample_groups) {
  # create mock Biobase::ExpressionSet
  eset_input = Biobase::ExpressionSet(
    assayData = x,
    featureData = Biobase::annotatedDataFrameFrom(x, byrow = T),
    protocolData = Biobase::annotatedDataFrameFrom(x, byrow = F)
  )
  rn = paste0("s", 1:nrow(x))
  Biobase::fData(eset_input) = data.frame(sequence = rn, row.names = rn, stringsAsFactors = T)
  Biobase::pData(eset_input) = data.frame(samples = colnames(x), condition = mask_sample_groups, row.names = colnames(x), stringsAsFactors = F)

  # actual normalization
  eset_norm = msEmpiRe::normalize(eset_input)
  x_norm = Biobase::exprs(eset_norm)
  x_norm[!is.finite(x_norm) | x_norm == 0] = NA
  return(x_norm)
}



#' sourced from MSqRob package, original code @ https://github.com/statOmics/MSqRob/blob/MSqRob0.7.6/R/preprocess_MaxQuant.R
#' @param exprs todo
#' @param weights todo
#'
#' @importFrom MASS rlm
normalize_rlr_MSqRob_implementation = function(exprs, weights = NULL) {
  mediandata = apply(exprs, 1, "median", na.rm = TRUE)
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



#' placeholder title
#' @param x todo
normalize_rlr = function(x) {
  y = normalize_rlr_MSqRob_implementation(x)
  dimnames(y) = dimnames(x)
  return(y)
}
