
#' Normalize a numerical matrix by the Variation Within, Mode Between (VWMB) algorithm
#'
#' @description
#' The normalization algorithm consists of two consecutive steps:
#' 1) samples are scaled within each group to minimize the overall `metric_within` among replicates
#' 2) summarize all samples per group by respective row mean values (from `row*sample` to a `row*group` matrix). Then rescale at the sample-group-level to minimize the overall `metric_between`
#'
#' The metric "var" will minimize the median of all variation-per-row.
#' The metric "mode" will minimize the modes of log-foldchange distribution between pairwise samples or groups (when set as `metric_within` and `metric_between` respectfully).
#'
#' Default setting is Variation Within, Mode Between (VWMB); `metric_within="var"` and `metric_between="mode"` will normalize such that median of variation-per-row is minimized per sample group, and subsequentially the log-foldchange distributions between all pairwise groups are minimized.
#'
#' Alternatively, one can also apply mode normalization within-group (Mode Within, Mode Between. MWMB) by setting `metric_within="mode"`. If the dataset has (unknown) covariates and a sufficient number of replicates, this might be beneficial because covariate-specific effects are not averaged out as they might be with `metric_within="var"`.
#'
#' To ignore groups and simply apply mode-between to all samples (not recommended!); `normalize_vwmb(x, groups=NA, metric_within="mode")`.
#'
#' To ignore groups and simply normalize all samples by reducing overall variation (not recommended!); `normalize_vwmb(x, groups=NA, metric_within="var")`.
#'
#' note; if you want to treat replicate samples that are flagged as 'exclude' upstream differently while still including them in the data matrix, you could set parameter `groups=paste(samples$group, samples$exclude)` to put them in separate groups.
#'
#' @examples \dontrun{
#' # Define a custom normalization function that we'll use in the MS-DAP pipeline later.
#' # All it does is apply the first part of the VWMB algorithm to reduce varation between replicates.
#' norm_within_only = function(x_as_log2, group_by_cols, ...) {
#'   normalize_vwmb(x_as_log2, groups = group_by_cols,
#'                  metric_within="var", metric_between = "")
#' }
#'
#' # now we use the main MS-DAP function as typical, but note that we here set the
#' # normalization algorithm to the function we defined in the lines above.
#' dataset = analysis_quickstart(
#'   dataset,
#'   norm_algorithm = "norm_within_only",
#'   # <other parameters here>
#' )
#' }
#' @param x numerical data matrix to normalize, should be log transformed
#' @param groups array describing the grouping of the columns in x (sample groups). Or alternatively set to NA to indicate there are no groups
#' @param metric_within how should replicate samples within a group be normalized? valid arguments: "var" reduce overall variation (default). "mode" reduce overall foldchange mode. pass empty string to disable
#' @param metric_between analogous to the metric_within parameter, how to normalize between groups? allowed parameters are "var" and "mode" (default). To disable, set groups=NA
#' @param include_attributes optionally, return some additional metrics as attributes of x. The "scaling" attribute describes the increase/decrease of each sample
#' @return normalized matrix x
#' @importFrom matrixStats colMedians rowSums2 rowMeans2
#' @export
normalize_vwmb = function(x, groups=NA, metric_within="var", metric_between="mode", include_attributes = FALSE) {
  #### input validation
  ## check parameters
  stopifnot(is.matrix(x) && typeof(x) == "double" && ncol(x) > 1 && nrow(x) > 1) # valid numeric matrix
  stopifnot((length(groups)==1 && is.na(groups)) || (length(groups) == ncol(x) && all(!is.na(groups)))) # groups should be NA or an array describing all columns of x
  stopifnot(length(metric_within)==1 && metric_within %in% c("mode","var", "")) # how to normalize within groups, with "" to disable
  stopifnot(length(metric_between)==1 && metric_between %in% c("mode","var", "")) # how to normalize between groups

  ## restrict group names: assume there are no groups if groups=NA  +  enforces character type for groups
  if(length(groups) == 1) {
    groups = rep("g1", ncol(x))
  } else {
    groups = paste0("g_", groups)
  }
  ugroups = unique(groups)
  group_count = length(ugroups)

  ## restrict input matrix
  # remove infinite values (eg; some upstream log(0)), then compute summary stats for input data
  x[!is.finite(x)] = NA
  x_median = stats::median(x, na.rm = T)
  if(include_attributes) {
    x_colmedians = matrixStats::colMedians(x, na.rm = T)
  }
  # replace colnames. if input matrix has colnames numerics-as-strings (eg; c(1,3,4)), indexing may be screwed downstream
  x_colnames = colnames(x)
  colnames(x) = paste0("V", 1:ncol(x))


  #### normalize within groups
  if(group_count > 1) {
    group_level_abundance = matrix(NA, nrow(x), group_count)
  }

  for (group_index in seq_along(ugroups)) { # group_index=1
    cols = which(groups == ugroups[group_index])

    # if there is just 1 sample in this group, no normalization to be done
    if (length(cols) == 1) {
      if(group_count > 1) {
        group_level_abundance[, group_index] = x[, cols]
      }
    } else {
      ## 1) feature selection
      rows_valid = matrixStats::rowSums2(!is.na(x[, cols])) >= 2 ## default. eg; rows with 1 feature are disregarded 'within group' anyway (no foldchange nor variation), so setting this filter also excludes them from between-group (for all groups with size >1)

      ## 2) scaling factor for each of 'cols'
      if(sum(rows_valid) >= 10) {
        s = NA
        if(metric_within == "mode") {
          s = norm_scales_fcmode(x[rows_valid, cols])
        }
        if(metric_within == "var") {
          s = norm_scales_var(x[rows_valid, cols])
        }

        ## 3) normalize all columns in this group
        if(all(!is.na(s))) {
          x[,cols] = x[,cols] + rep(s, rep.int(nrow(x), length(cols))) # slightly faster than looping columns
          # for(index in seq_along(cols)) { x[, cols[index]] = x[ , cols[index]] + s[index] }
        }
      } else {
        # this should never occur for any reasonably sized dataset
        warning(paste("need at least 10 rows with a value, not normalizing these samples;", paste(x_colnames[cols], collapse=", ")))
      }

      ## 4) store group-level values for downstream between-group normalization
      # mean value for each row represents the group-level abundance when scaling between groups
      if(group_count > 1) {
        group_level_abundance[rows_valid, group_index] = matrixStats::rowMeans2(x[rows_valid, cols], na.rm=T)
      }
    }
  }


  #### normalize between groups
  if(group_count > 1 && metric_between != "") {
    # analogous to above
    # metric_between is either "mode" or "var"
    if(metric_between == "var") {
      s = norm_scales_var(group_level_abundance)
    } else {
      s = norm_scales_fcmode(group_level_abundance)
    }
    for (group_index in seq_along(ugroups)) {
      cols = which(groups == ugroups[group_index])
      x[, cols] = x[, cols] + s[group_index]
    }
  }


  #### finally, scale entire result matrix such that the median intensity is the same as input (deals with the rare case where abundances in one group are much higher, and this normalization consequentially shifts overall abundance levels @ between-group scaling)
  x = x + (x_median - stats::median(x, na.rm = T))
  # reset colnames
  colnames(x) = x_colnames
  # threshold very low values
  x = threshold_numerics(x, 1)


  #### optionally, add some metadata for downstream QC
  if(include_attributes) {
    # scaling factors per column
    x_scaling = x_colmedians - matrixStats::colMedians(x, na.rm = T)
    x_scaling[!is.finite(x_scaling)] = 0
    attr(x, "scaling") = x_scaling
    # features used to normalize between groups
    if(group_count > 1) {
      attr(x, "features_between") = group_level_abundance
    }
  }

  return(x)
}



#' Given some peptide-level matrix, normalize samples using foldchange mode-between at protein-level
#'
#' @param x numerical data matrix with log2 peptide abundances
#' @param groups array describing the grouping of the columns in x
#' @param proteins array describing the grouping of the rows in x
#' @param rollup_algorithm algorithm for combining peptides to protein-level data, on which in subsequent steps the normalization/scaling factors will be computed. Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#' @return normalized matrix x
#' @export
normalize_modebetween_protein = function(x, proteins, groups, rollup_algorithm = "maxlfq") {
  ## input validation
  stopifnot(is.matrix(x) && typeof(x) == "double" && ncol(x) > 1 && nrow(x) > 1)
  stopifnot(length(proteins) == nrow(x) && all(!is.na(proteins)) && (all(is.character(proteins) & proteins != "") || all(is.integer(proteins))) )
  stopifnot(length(groups) == ncol(x) && all(!is.na(groups)) && (all(is.character(groups) & groups != "") || all(is.integer(groups))) )
  # rollup_algorithm parameter is validated in rollup function

  ## rollup
  # from peptide matrix to long-format tibble
  tib = as_tibble(x) %>%
    add_column(peptide_id = 1:nrow(x),
               protein_id = proteins) %>%
    tidyr::pivot_longer(cols = c(-protein_id, -peptide_id), names_to = "sample_id", values_to = "intensity") %>%
    filter(is.finite(intensity) & intensity > 0)
  # rollup
  mat_protein = rollup_pep2prot(tib, intensity_is_log2 = TRUE, rollup_algorithm = rollup_algorithm, return_as_matrix = TRUE)

  ## importantly, enforce consistent column order to preserve alignment with sample-to-group assignments
  stopifnot(ncol(x) == ncol(mat_protein) && all(colnames(x) %in% colnames(mat_protein)))
  # use match() instead of direct key/string-based indexing because some samples may have names like 1,2,3,4
  # (eg; if key_sample is used for column names instead of sample_id, as we do in filter_dataset() )
  mat_protein = mat_protein[,match(colnames(x), colnames(mat_protein)),drop=F]

  ## apply VWMB function to obtain mode-between scaling factors
  mat_protein = normalize_vwmb(mat_protein, groups = groups, metric_within = "", metric_between = "mode", include_attributes = TRUE) # disable within-group normalization
  scale_per_sample = attr(mat_protein, "scaling")

  # apply scaling factors, determined through normalization on protein-level, to peptide-level (simple loop is more efficient than transpose+scale+transpose)
  for(j in seq_along(scale_per_sample)) {
    x[,j] = x[,j] - scale_per_sample[j]
  }
  return(x)
}



#' compute the mode of all sample/column pairs
#' column pairs that have less than 10 values in common are skipped (result value left at default 0)
#' @param x log-transformed data matrix, where columns are samples and rows are features
#' @param min_overlap minimum number of rows that have a finite foldchange for any given combination of columns (otherwise, mode foldchange is not computed)
#' @return a matrix of ncol*ncol pairwise foldchange-modes
pairwise_modes = function(x, min_overlap = 10) {
  m = matrix(0, ncol(x), ncol(x), dimnames = list(colnames(x), colnames(x)))
  for (i in 1:ncol(x)) { # i=1;j=2
    for (j in 1:ncol(x)) {
      if (i > j) {
        # foldchanges from column i to j
        fc = x[, i] - x[, j]
        fc = fc[is.finite(fc)]
        # don't do anything if less than N data points
        if(length(fc) >= min_overlap) {
          # find the mode of all (finite) log foldchanges
          m[i, j] = get_mode(fc)
          # from j to i is the opposite
          m[j, i] = -1 * m[i, j]
        }
      }
    }
  }
  return(m)
}



#' find the mode in a numeric array using the density function from base R
#' @param x numeric array
#' @return mode of x
get_mode = function(x) {
  if(length(x) == 0) return(NA)
  if(length(x) == 1) return(x)

  density_estimate = density(x, na.rm = T) # default density kernel (left as-is for consistency with prior MS-DAP versions)
  # density_estimate = density(x, bw = "sj", adjust = 0.9, na.rm = T) # alternative bandwidth; we need to benchmark to validate this is an improvement for the vast majority of datasets (prior to updating previous default)
  # which.max() returns first 'max' value (thus always 1 return value). eg; which.max(c(1,1,2,2,1))
  return(density_estimate$x[which.max(density_estimate$y)])
}



#' scale sample k by scaling factor s_k  -->>  don't have to re-compute all scalings, can just update the pairwise matrix
#' @param m matrix with pairwise foldchange modes
#' @param s array of scaling factors for each sample
#' @return updated matrix m
adjust_modes = function(m, s) {
  for(k in seq_along(s)) {
    m[k,] = m[k,] + s[k]
    m[,k] = m[,k] - s[k]
  }
  return(m)
}



#' minimize the foldchange-mode between all columns
#' @param x log-transformed numeric matrix
#' @param min_overlap minimum number of rows that have a finite foldchange for any given combination of columns (otherwise, mode foldchange is not computed)
#' @return scaling factor for each column in x
#' @importFrom matrixStats rowSums2
norm_scales_fcmode = function(x, min_overlap = 10) {
  ### compute the foldchange-mode between each pair of columns
  fcm = pairwise_modes(x, min_overlap = min_overlap)

  ### apply MLE to scale all samples such that overall, the foldchange-mode between any pair of samples is minimal
  # helper function for MLE
  mle_norm = function(...) {
    m = adjust_modes(fcm, unlist(...))
    return( sum(abs(m)) )
  }

  # apply MLE with default settings, starting point is normalization to some reference sample
  # reference sample = least overall distance to other samples
  row_most_similar = order(matrixStats::rowSums2(abs(fcm), na.rm = T), decreasing = F, na.last = T)[1]
  s_mle = optim(par = fcm[row_most_similar,], fn = mle_norm, method = "BFGS")

  # TODO: optionally, we could print some some metrics on the chosen scaling (how many MLE iterations, overall score before/after)
  return(s_mle$par)
}



#' minimize variation over all rows in a numeric matrix
#' @param x log-transformed numeric matrix
#' @return scaling factor for each column in x
#' @importFrom matrixStats rowMeans2 rowSds colMedians
norm_scales_var = function(x) {
  x_scaled = x - matrixStats::rowMeans2(x, na.rm = T)
  # apply scaling to each sample, then compute standard deviation on each row; y = t(t(x_scaled) + s)
  mle_norm = function(...) {
    s = unlist(...)
    # use the median instead of the mean, as the latter is too sensitive to outliers (eg; upstream data has some errors so few features have an extreme value -->> dominates the scaling of all features)
    return( stats::median(matrixStats::rowSds(x_scaled + rep(s, rep.int(nrow(x_scaled), ncol(x_scaled))), na.rm=T), na.rm=T) )
    # optimizations; column-wise matrix operations; https://stackoverflow.com/a/32364355   https://stats.stackexchange.com/a/51750
    # optimizations; median @ https://stackoverflow.com/questions/34771088/why-is-standard-r-median-function-so-much-slower-than-a-simple-c-alternative
  }
  # apply MLE with default settings, starting point is scaling such that median values are the same
  x_scaled_colmedian = matrixStats::colMedians(x_scaled, na.rm=T)
  s_init = mean(x_scaled_colmedian) - x_scaled_colmedian
  s_mle = optim(par = s_init, fn = mle_norm, method = "BFGS")

  # TODO: optionally, we could print some some metrics on the chosen scaling (how many MLE iterations, overall score before/after)
  return(s_mle$par)
}



### rudimentatry testing: validate that we properly deal with a lack of overlap between samples or groups
# xx = rnorm(20)
# x = cbind(c(xx, rep(NA, length(xx))),
#           c(xx+1, rep(NA, length(xx))),
#           c(rep(NA, length(xx)), xx+3),
#           c(rep(NA, length(xx)), xx+4))
#
# # plot data as-is
# pheatmap::pheatmap(x, na_col = "grey", cluster_cols = F, cluster_rows = F)
#
# # normalize within replicates and between groups (no overlap between groups)
# y = normalize_vwmb(x, groups = c(1,1,2,2))
# pheatmap::pheatmap(y, na_col = "grey", cluster_cols = F, cluster_rows = F)
#
# # normalize by 'reducing overall variation', while some samples have zero overlap
# y = normalize_vwmb(x, groups = NA)
# pheatmap::pheatmap(y, na_col = "grey", cluster_cols = F, cluster_rows = F)
#
# # define groups such that there is zero overlap between replicates. results in a warning
# y = normalize_vwmb(x, groups = c(1,2,1,2))
# pheatmap::pheatmap(y, na_col = "grey", cluster_cols = F, cluster_rows = F)
###
