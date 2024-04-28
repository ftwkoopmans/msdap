
#' Compare the number of peptide detection counts between groups
#'
#' The computation of scores is detailed in the online vignitte for "differential testing".
#'
#' @param dataset a valid dataset
#' @param min_peptides_observed minimum number of peptides for a protein to pass filtering rules (i.e. otherwise, no z-score is computed)
#' @param min_samples_observed minimum number of samples where a protein should be observed with at least `min_peptides_observed` peptides (in either group) when comparing a contrast of group A vs B
#' @param min_fraction_observed for differential detection only; analogous to `diffdetect_min_samples_observed`, but here you can specify the fraction of samples where a protein needs to be detected in either group (within the respective contrast). default; 0.5 (50% of samples)
#' @param count_mode whether z-scores for only detected peptides, or all observed peptides, should be computed. options; "auto", "detect", "quant"
#' @param rescale_counts_per_sample boolean, indicating whether per sample, detect counts should be rescaled (i.e. to normalize across samples for the number of detected peptides)
#' @param return_wide_format boolean, return table in simplified wide-format (TRUE) or return all data in long-format table (FALSE, default)
#' @export
differential_detect = function(dataset, min_peptides_observed = 1L, min_samples_observed = 3, min_fraction_observed = 0.5, count_mode = "auto", rescale_counts_per_sample = TRUE, return_wide_format = FALSE) {

  # helper function
  summarize_proteins = function(mat, mat_scaled, sid1, sid2, min_peptides_observed, min_samples_observed, min_fraction_observed, rescale_counts_per_sample) {
    mat1 = mat[,colnames(mat) %in% sid1, drop=F]
    mat2 = mat[,colnames(mat) %in% sid2, drop=F]
    ### proteins with valid data for computing a z-score;
    ### rows that were detected with N+ peptides in M+ samples  (in either experimental condition @ current contrast)
    # count number of columns/samples with at least N peptides
    mat1_rowcount = matrixStats::rowSums2(mat1 >= min_peptides_observed)
    mat2_rowcount = matrixStats::rowSums2(mat2 >= min_peptides_observed)
    # criteria per sample group / experimental condition; hardcoded minimum number of samples & fraction-of-samples
    threshold1 = max(min_samples_observed, ceiling(ncol(mat1) * min_fraction_observed))
    threshold2 = max(min_samples_observed, ceiling(ncol(mat2) * min_fraction_observed))
    rows_pass = mat1_rowcount >= threshold1 | mat2_rowcount >= threshold2

    ### collect data for z-score and output table
    # for each condition; nobs, nobs_scaled, upep_max, usample_max
    z = tibble::tibble(
      protein_id = rownames(count_detect),
      nsample1 = matrixStats::rowSums2((mat1 != 0L) + 0L), # from bool to int prior to matrixStats
      nsample2 = matrixStats::rowSums2((mat2 != 0L) + 0L),
      npep1 = matrixStats::rowMaxs(mat1),
      npep2 = matrixStats::rowMaxs(mat2),
      npep_max = matrixStats::rowMaxs(mat),
      nobs1 = matrixStats::rowSums2(mat1),
      nobs2 = matrixStats::rowSums2(mat2),
      nobs1_scaled = matrixStats::rowSums2(mat_scaled[,colnames(mat_scaled) %in% sid1, drop=F]),
      nobs2_scaled = matrixStats::rowSums2(mat_scaled[,colnames(mat_scaled) %in% sid2, drop=F]),
      pass_filters = rows_pass,
      log2fc = NA,
      zscore = NA,
      fracdiff = NA
    )

    # user doesn't want to use the normalized counts. we can simply overwrite
    if(rescale_counts_per_sample == FALSE) {
      z$nobs1_scaled = z$nobs1
      z$nobs2_scaled = z$nobs2
    }

    ### compute z-score for all proteins that passed, but retain all other proteins in result table (more flexible for downstream use of these count data)
    if(any(z$pass_filters)) {
      # log2 scaled foldchange, padding each condition with minimum non-zero score (to deal with zero cases).
      # in most datasets, this'll be virtually the same as 1
      min1 = min2 = 1
      rows1 = z$nobs1_scaled > 0 & z$pass_filters
      rows2 = z$nobs2_scaled > 0 & z$pass_filters
      if(any(rows1)) min1 = min(z$nobs1_scaled[rows1])
      if(any(rows2)) min2 = min(z$nobs2_scaled[rows2])

      if(any(rows1) || any(rows2)) {
        ### log2FC = log2(b+1) - log2(a+1)   (but since we rescaled the count matrix, we use the minimum observed value instead of hardcoded +1. should not make much difference)
        z$log2fc = suppressWarnings( log2(z$nobs2_scaled + min2) - log2(z$nobs1_scaled + min1) )
        z$log2fc[!is.finite(z$log2fc) | z$pass_filters == FALSE] = NA # enforce NA's over inf/NaN
        ### standardize score
        ## updated metric in MS-DAP release 1.6
        # MS-DAP release 1.6.1  -->>  center log2fc in results as well  (previously we centered log2fc prior to z-score computation, but didn't return centered data)
        z$log2fc = z$log2fc - stats::median(z$log2fc, na.rm=T)
        sd_rescale = 1 # init variable, overwritten in all cases
        fit_sd = suppressWarnings(fit_t_dist_fixed_mu(z$log2fc))
        if(length(fit_sd) == 2 && is.finite(fit_sd[1])) {
          sd_rescale = fit_sd[1]
        } else {
          sd_rescale = sd(z$log2fc, na.rm = TRUE)
        }
        z$zscore = suppressWarnings(z$log2fc / sd_rescale)
        # z$zscore = suppressWarnings((z$log2fc - get_mode(z$log2fc)) / sd_rescale)
        # z$zscore = suppressWarnings((z$log2fc - get_mode(z$log2fc)) / sd(z$log2fc, na.rm = T))
        # z$zscore = suppressWarnings((z$log2fc - get_mode(z$log2fc)) / mad(z$log2fc, na.rm = T))
        ## previous metric; plain z-score  (but vulnerable to distributions that actually have relatively many real diff values, e.g. IP of WT+KO where 25/125 proteins are actually different = wide t-distribution)
        # z$zscore = suppressWarnings((z$log2fc - mean(z$log2fc)) / sd(z$log2fc, na.rm = T))
        z$zscore[!is.finite(z$zscore)] = NA
        ### "proportion of change" metrix as a proxy for absolute diff that accounts for number of peptides; (absolute diff / number of peptides) / number of samples
        # ! uses the raw counts, not the 'normalized' counts used for z-scores
        z$fracdiff = suppressWarnings( ((z$nobs2 - z$nobs1) / z$npep_max)  /  max(ncol(mat1), ncol(mat2)) )
        z$fracdiff[!is.finite(z$fracdiff)] = NA
      }
    }

    return(z)
  }



  # remove preexisting results
  dataset$dd_proteins = NULL

  ### input validation
  if(length(count_mode) != 1 || ! count_mode %in% c("auto", "detect", "quant") ) {
    append_log("differential detection parameter 'count_mode' must be a single string, with valid options; 'auto', 'detect', 'quant'", type = "error")
  }
  if(length(min_peptides_observed) != 1 || (!is.na(min_peptides_observed) && !(is.numeric(min_peptides_observed) && is.finite(min_peptides_observed)) ) ) {
    append_log("differential detection parameter 'min_peptides_observed' must be a single integer value (or NA to skip this analysis)", type = "error")
  }
  if(length(min_samples_observed) != 1 || (!is.na(min_samples_observed) && !(is.numeric(min_samples_observed) && is.finite(min_samples_observed)) ) ) {
    append_log("differential detection parameter 'min_samples_observed' must be a single integer value (or NA to skip this analysis)", type = "error")
  }
  if(length(min_fraction_observed) != 1 || !is.finite(min_fraction_observed) || min_fraction_observed < 0 || min_fraction_observed > 1) {
    append_log("differential detection parameter 'min_fraction_observed' must be a single numeric value between 0 and 1", type = "error")
  }

  # user doesn't want to perform differential detect; min_peptides_observed or min_samples_observed is set to NA
  if(is.na(min_peptides_observed) || is.na(min_samples_observed)) {
    append_log("differential detection parameters set to NA, this analysis is cancelled", type = "warning")
    return(dataset)
  }

  if(!is.na(min_peptides_observed) && !is.integer(min_peptides_observed)) {
    min_peptides_observed = as.integer(ceiling(min_peptides_observed))
  }
  if(!is.na(min_samples_observed) && !is.integer(min_samples_observed)) {
    min_samples_observed = as.integer(ceiling(min_samples_observed))
  }

  if("fraction" %in% colnames(dataset$samples) && n_distinct(dataset$samples$fraction) > 1) {
    append_log("you have to merge fractionated samples prior to differential detection. Typically, you want to simply use the analysis_quickstart() to process the dataset first. Advanced users may call merge_fractionated_samples(dataset) manually.", type = "error")
  }

  column_contrasts = dataset_contrasts(dataset)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential detection analysis is cancelled", type = "warning")
    return(dataset)
  }


  append_log(sprintf("differential detection analysis: min_samples_observed=%d min_fraction_observed=%.2f", min_samples_observed, min_fraction_observed), type = "info")
  only_detected_peptides = FALSE
  if(count_mode == "auto") {
    # for DIA data, count only peptides that are confidently detected (i.e. disregard MBR)
    only_detected_peptides = is_dia_dataset(dataset)
  } else {
    only_detected_peptides = (count_mode == "detect")
  }



  ### protein count data from the dataset$peptides tibble

  # subset of the peptide table where we have intensity values (i.e. not decoy)  AND  ignore 'exclude' samples
  x = dataset$peptides %>% filter(is.finite(intensity) & sample_id %in% (dataset$samples %>% filter(exclude == F) %>% pull(sample_id)) )

  # peptide*sample detect matrix
  tmp = x %>% select(peptide_id, sample_id, detect) %>% tidyr::pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "detect")
  rm(x)
  m = as.matrix(tmp %>% select(-peptide_id))
  rownames(m) = tmp$peptide_id
  rm(tmp) # m[1:3, 1:4]
  m__protein_id = dataset$peptides$protein_id[data.table::chmatch(rownames(m), dataset$peptides$peptide_id)]
  # detected peptide are TRUE, quant but not detected are FALSE (e.g. MBR), values not observed in a sample are NA
  # here we convert to integer scores 2, 1, 0, respectively
  m = m + 1
  m[is.na(m)] = 0
  mode(m) = "integer"
  m_detect = m == 2L
  m_quant = m >= 1L
  rm(m)
  # note that for some datasets, only "detected" peptides are present so 'detect' and 'quant' counts are the same
  compute_quant_counts = only_detected_peptides == FALSE && any(m_detect != m_quant)

  # normalize
  sample_count_detect = matrixStats::colSums2(m_detect)
  sample_count_detect__scaling_factor = sample_count_detect / mean(sample_count_detect)

  # protein*sample count matrices
  df_detect = stats::aggregate(m_detect, by = list(protein_id = m__protein_id), FUN = sum) # slow !
  count_detect = as.matrix(df_detect[,-1])
  rownames(count_detect) = df_detect$protein_id
  rm(df_detect)
  # rescale by total count events per sample
  count_detect__scaled = count_detect
  for(j in 1L:ncol(count_detect)) {
    count_detect__scaled[,j] = count_detect[,j] / sample_count_detect__scaling_factor[j]
  }

  # analogous
  if(compute_quant_counts) {
    sample_count_quant = matrixStats::colSums2(m_quant)
    sample_count_quant__scaling_factor = sample_count_quant / mean(sample_count_quant)
    df_quant = stats::aggregate(m_quant, by = list(protein_id = m__protein_id), FUN = sum)
    count_quant = as.matrix(df_quant[,-1])
    rownames(count_quant) = df_quant$protein_id
    rm(df_quant)
    count_quant__scaled = count_quant
    for(j in 1L:ncol(count_quant)) {
      count_quant__scaled[,j] = count_quant[,j] / sample_count_quant__scaling_factor[j]
    }
  }



  ### compare counts within each contrast

  tib_results = NULL
  for(col_contr in column_contrasts) { # col_contr = column_contrasts[1]
    ### sample_id per condition ('side of each contrast')
    contr_samples = dataset$samples %>%
      select(sample_id, exclude, condition=!!col_contr) %>%
      filter(exclude == FALSE & condition != 0) %>%
      arrange(condition)
    # set of sample_id for condition 1 and 2
    sid1 = contr_samples %>% filter(condition == 1) %>% pull(sample_id)
    sid2 = contr_samples %>% filter(condition == 2) %>% pull(sample_id)

    # compute summary stats and z-score for each protein_id
    z_detect = summarize_proteins(count_detect, count_detect__scaled, sid1, sid2, min_peptides_observed, min_samples_observed, min_fraction_observed, rescale_counts_per_sample)
    # append to results
    tib_results = bind_rows(tib_results, z_detect %>% mutate(contrast = col_contr, type = "detect"))

    # analogous
    if(compute_quant_counts) {
      z_quant = summarize_proteins(count_quant, count_quant__scaled, sid1, sid2, min_peptides_observed, min_samples_observed, min_fraction_observed, rescale_counts_per_sample)
      tib_results = bind_rows(tib_results, z_quant %>% mutate(contrast = col_contr, type = "quant"))
    }
  }

  if(is.data.frame(tib_results)) {
    tib_results = tib_results %>% mutate_all(unname)
  }

  # optionally, user can return results as a matrix
  if(return_wide_format == TRUE && nrow(tib_results) > 0) {
    tmp = tib_results # %>% filter(!is.na(zscore))
    if(nrow(tmp) > 0) {
      dataset$dd_proteins = tmp %>%
        select(protein_id, contrast, type, npep = npep_max, nobs1, nobs2, log2fc, fracdiff, zscore) %>%
        tidyr::pivot_wider(id_cols = c("protein_id", "contrast"), names_from = "type", names_prefix = "count_", values_from = c("npep", "nobs1", "nobs2", "log2fc", "fracdiff", "zscore"))
    }
  } else {
    dataset$dd_proteins = tib_results
  }

  return(dataset)
}



#' convert the results from differential detection analysis from a long-format tibble to wide-format
#'
#' @param dataset your dataset. if 'dd_proteins' is lacking, result is empty tibble
#' @importFrom tidyr pivot_wider everything
diffdetect_results_to_wide = function(dataset) {
  if(!is_tibble(dataset$dd_proteins) || nrow(dataset$dd_proteins) == 0) {
    return(tibble())
  }

  ## collect number of 'detect events' per sample group (outlier samples not taken into account)
  # bugfix
  # account for default behaviour by tidyr::pivot_wider()  -->>  if there is just 1 "values_from" column, the names are not prefixed but if there are multiple (eg. both detect and quant scores), the input column names are prefixed
  # we always want this "prefix" because in current implementation, for DIA we only have "peptide detect count" based z-scores while for DDA we get both detect and quant z-scores
  # solution: `names_glue` parameter

  tib_counts = dataset$peptides %>%
    select(protein_id, sample_id, detect) %>%
    left_join(dataset$samples %>% select(sample_id, group), by = "sample_id") %>%
    group_by(protein_id, group) %>%
    summarise(count_peptides_detected_within_group = sum(detect), count_peptides_quantified_within_group = n()) %>%
    ungroup() %>%
    tidyr::pivot_wider(id_cols = "protein_id", names_from = "group", names_glue = "{.value}__{group}",
                       values_from = c("count_peptides_detected_within_group", "count_peptides_quantified_within_group"),
                       values_fill = list(count_peptides_detected_within_group = 0, count_peptides_quantified_within_group = 0) )


  tib_zscores = dataset$dd_proteins %>%
    select(protein_id, contrast, type, zscore) %>%
    filter(is.finite(zscore)) %>%
    tidyr::pivot_wider(id_cols = "protein_id", names_from = c("type", "contrast"), names_glue = "{.value}_{type}__{contrast}", values_from = "zscore")


  return( full_join(tib_counts, tib_zscores, by = "protein_id") )
}
