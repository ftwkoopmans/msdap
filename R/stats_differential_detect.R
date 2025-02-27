
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
  start_time = Sys.time()

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

  if(!is.list(dataset$contrasts) || length(dataset$contrasts) == 0) {
    append_log("no contrasts have been defined, differential detection analysis is cancelled", type = "warning")
    return(dataset)
  }

  append_log(sprintf("differential detection analysis: min_peptides_observed=%d min_samples_observed=%d min_fraction_observed=%.2f count_mode=%s rescale_counts_per_sample=%s",
                     min_peptides_observed, min_samples_observed, min_fraction_observed, count_mode, rescale_counts_per_sample), type = "info")
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
  sample_count_detect = matrixStats::colSums2(m_detect) # number of detects per sample/column
  tmp = sample_count_detect[sample_count_detect > 0] # non-zero values
  tmp2 = ifelse(length(tmp) == 0, 1, mean(tmp)) # mean of non-zero values, or 1 if there are none
  sample_count_detect__scaling_factor = sample_count_detect / tmp2
  rm(tmp, tmp2)
  # if a sample has 0 detected peptides, the scaling factor will also be 0 which will cause NA's downstream
  sample_count_detect__scaling_factor[!is.finite(sample_count_detect__scaling_factor) | sample_count_detect__scaling_factor <= 0.001] = 1 # don't scale samples with invalid scaling factors

  # protein*sample count matrices
  count_detect = matrix_grouped_column_aggregate(m_detect, m__protein_id, FUN = sum)
  # rescale by total count events per sample
  count_detect__scaled = count_detect
  for(j in 1L:ncol(count_detect)) {
    count_detect__scaled[,j] = count_detect[,j] / sample_count_detect__scaling_factor[j]
  }

  # analogous
  if(compute_quant_counts) {
    sample_count_quant = matrixStats::colSums2(m_quant)
    tmp = sample_count_quant[sample_count_quant > 0] # non-zero values
    tmp2 = ifelse(length(tmp) == 0, 1, mean(tmp)) # mean of non-zero values, or 1 if there are none
    sample_count_quant__scaling_factor = sample_count_quant / tmp2
    rm(tmp, tmp2)
    # if a sample has 0 detected peptides, the scaling factor will also be 0 which will cause NA's downstream
    sample_count_quant__scaling_factor[!is.finite(sample_count_quant__scaling_factor) | sample_count_quant__scaling_factor <= 0.001] = 1 # don't scale samples with invalid scaling factors

    count_quant = matrix_grouped_column_aggregate(m_quant, m__protein_id, FUN = sum)
    count_quant__scaled = count_quant
    for(j in 1L:ncol(count_quant)) {
      count_quant__scaled[,j] = count_quant[,j] / sample_count_quant__scaling_factor[j]
    }
  }



  ### compare counts within each contrast

  tib_results = NULL
  for(contr in dataset$contrasts) {
    # set of sample_id for condition 1 and 2
    sid1 = contr$sampleid_condition1
    sid2 = contr$sampleid_condition2

    # compute summary stats and z-score for each protein_id
    z_detect = summarize_proteins(count_detect, count_detect__scaled, sid1, sid2, min_peptides_observed, min_samples_observed, min_fraction_observed, rescale_counts_per_sample)
    # append to results
    tib_results = bind_rows(tib_results, z_detect %>% mutate(contrast = contr$label, type = "detect"))

    # analogous
    if(compute_quant_counts) {
      z_quant = summarize_proteins(count_quant, count_quant__scaled, sid1, sid2, min_peptides_observed, min_samples_observed, min_fraction_observed, rescale_counts_per_sample)
      tib_results = bind_rows(tib_results, z_quant %>% mutate(contrast = contr$label, type = "quant"))
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

  append_log_timestamp("differential detection analysis", start_time)
  return(dataset)
}



#' Protein differential detection filtering using a hardcoded number of peptide observations per condition
#'
#' @description
#' ## Pseudocode
#'
#' Number of peptides differentially detected in condition 1 as compared to condition 2:
#'
#' `npep_pass1 = sum(n1 - n2 >= k1_diff)`
#'
#' Where n1 and n2 are vectors with the number of detects per peptide in conditions 1 and 2, respectively,
#' and k1_diff is the sample count threshold for differential detection.
#'
#' Now we can test at protein-level, directionally!
#' We're looking for proteins that are near-exclusively detected in condition 1 and their nobs ratio is in the same direction.
#'
#' `npep_pass1 >= npep_pass & nobs1 > nobs_ratio * nobs2`
#'
#' Where npep_pass and nobs_ratio are user-provided thresholds, nobs1 is the sum of n1 (see further return value specification for this function).
#' Above rule tests enrichment in condition 1, afterwards we also test analogously with conditions 1 and 2 swapped.
#'
#'
#' ## real-world example data:
#'
#' 2 conditions, 6 replicates each:
#'
#' \tabular{llrr}{
#' \strong{protein_id} \tab \strong{peptide_id} \tab \strong{n1} \tab \strong{n2} \cr
#' P09041 \tab ALMDEVVK_2 \tab 2 \tab 6 \cr
#' P09041 \tab LTLDKVDLK_2 \tab 1 \tab 5
#' }
#'
#' n2-n1 diff is 4 for both peptides, nobs ratio was 3.7-fold enrichment in condition 2.
#'
#' @param dataset your dataset. At least 1 contrast should have been specified prior
#' @param k_diff peptide-level criterium: peptide must be detected in at least k more samples in condition 1 versus condition 2, or vice versa
#' @param frac_diff peptide-level criterium: peptide must be detected in at least x% more samples in condition 1 versus condition 2, or vice versa
#' @param npep_pass minimum number of peptides that must pass the peptide-level differential detection criteria. Default: 2
#' @param nobs_ratio minimum enrichment ratio between experimental conditions for the total peptide detection count per protein (see output table specification, `nobs1` and `log2fc_nobs2/nobs1` to which this filter is applied). Note that this value is NOT on log2 scale, i.e. set 2 for 2-fold enrichment. Default: 4
#' @param int_ratio analogous to `nobs_ratio`, but for the enrichment in sum peptide intensity values. Default: 0 (disabled)
#' @param normalize_intensities normalize the protein intensity matrix prior to computing sum intensities and intensity ratios. Default: TRUE
#' @examples \dontrun{
#' ## example 1:
#' # default / stringent: find proteins that have at least 2 peptides
#' # detected in 66% more samples in either condition (with a minimum of 3 samples)
#' # AND the overall detection rate is larger than 3-fold
#' x = differential_detection_filter(
#'   dataset, k_diff = 3, frac_diff = 0.66, npep_pass = 2, nobs_ratio = 3
#' )
#'
#' # code snippet to create a pretty-print table with resulting proteins
#' y = x %>%
#'   # only retain proteins that match input criteria
#'   filter(pass) %>%
#'   # add protein metadata and rearrange columns
#'   left_join(dataset$proteins %>%
#'     select(protein_id, fasta_headers, gene_symbols_or_id),
#'     by = "protein_id") %>%
#'   select(contrast, protein_id, fasta_headers, gene_symbols_or_id,
#'          tidyselect::everything()) %>%
#'   # for prettyprint, trim the contrast names
#'   mutate(contrast = gsub(" *#.*", "", sub("^contrast: ", "", contrast))) %>%
#'   # sort data by column, then by ratio
#'   arrange(contrast, `log2fc_nobs2/nobs1`)
#' print(y)
#'
#' ## example 2:
#' # a more lenient filter: apply criteria to only 1 (or 0) peptides but
#' # rely mostly on the nobs_ratio criterium and additionally
#' # add post-hoc filtering on the total number of peptides per protein
#' # and require either protein to have an overall 50% detection rate
#' # (across peptides and samples)
#' x = differential_detection_filter(
#'   dataset, nobs_ratio = 3, npep_pass = 1, k_diff = 3, frac_diff = 0.66
#' ) %>%
#'   mutate(pass = pass & npep_total > 1 & pmax(fracobs1, fracobs1) >= 0.5)
#' }
#' @returns A tibble where each row describes 1 proteingroup ("protein_id") in 1 contrast, with the following columns:
#'
#' - `npep_total` = total number of peptides detected across any of the samples in the current contrast. Useful for post-hoc filtering, e.g. when you do not set stringent criteria for `npep_pass`
#' - `npep_pass1` = number of peptides that pass the specified filtering rules in condition 1 of the current contrast
#' - `npep_pass2` = analogous to `npep_pass1`, but for condition 2
#' - `nobs1` = sum of peptide detections across all samples in condition 1 (i.e. each detected peptide is counted once per sample, this is the total sum across peptides*samples for respective protein_id and condition). Useful for post-hoc filtering, e.g. when you do not set stringent criteria for `npep_pass`
#' - `nobs2` = analogous to `nobs1`, but for condition 2
#' - `fracobs1` = percentage of all possible detections made within condition 1 = number of observed datapoints / (#peptide * #sample)
#' - `fracobs2` = analogous to `fracobs1`, but for condition 2
#' - `log2fc_nobs2/nobs1` = log2 foldchange of observation counts. Positive values are enriched in condition 2. Proteins exclusive to condition 2 have value `Inf` and exclusive to condition 1 have value `-Inf`
#' - `int1` = sum peptide intensity across all samples in condition 1
#' - `int2` = sum peptide intensity across all samples in condition 2
#' - `log2fc_int2/int1` = log2 foldchange of protein intensities, analogous to `log2fc_nobs2/nobs1`
#' - `pass` = protein matches input filtering
#' @export
differential_detection_filter = function(dataset, k_diff = NA, frac_diff = NA, npep_pass = 2L, nobs_ratio = 3, int_ratio = 0, normalize_intensities = TRUE) {
  # input validation
  has_kdiff = length(k_diff) == 1 && is.numeric(k_diff) && !is.na(k_diff) && k_diff >= 0
  has_fracdiff = length(frac_diff) == 1 && is.numeric(frac_diff) && !is.na(frac_diff) && frac_diff >= 0 && frac_diff <= 1
  stopifnot("k_diff and frac_diff parameters are both invalid" = has_kdiff | has_fracdiff)
  stopifnot("nobs_ratio parameter is invalid" = length(nobs_ratio) == 1 && is.numeric(nobs_ratio) && !is.na(nobs_ratio) && nobs_ratio >= 0)
  stopifnot("int_ratio parameter is invalid" = length(int_ratio) == 1 && is.numeric(int_ratio) && !is.na(int_ratio) && int_ratio >= 0)
  stopifnot("npep_pass parameter is invalid" = length(npep_pass) == 1 && is.numeric(npep_pass) && !is.na(npep_pass) && npep_pass >= 0)
  if (!is.list(dataset$contrasts) || length(dataset$contrasts) == 0) {
    append_log("no contrasts have been defined, differential detection analysis is cancelled", type = "warning")
    return(NULL)
  }

  # round up and enforce integer types for some params
  npep_pass = as.integer(ceiling(npep_pass))
  if(has_kdiff) {
    k_diff = as.integer(ceiling(k_diff))
  }

  # prepare count table for faster computation downstream
  count_table = dataset$peptides %>%
    filter(detect == TRUE) %>%
    select(protein_id, peptide_id, sample_id, detect) %>%
    mutate(detect = as.integer(detect)) %>%
    pivot_wider(id_cols = c("protein_id", "peptide_id"), names_from = "sample_id", values_from = "detect", values_fill = 0L)

  # in matrix format for fast summarization
  count_matrix = as.matrix(count_table %>% select(-protein_id, -peptide_id))
  stopifnot(is.finite(count_matrix) & count_matrix %in% c(0L, 1L))


  # analogous; intensity table
  abundance_table = dataset$peptides %>%
    filter(is.finite(intensity)) %>%
    select(protein_id, peptide_id, sample_id, intensity) %>%
    pivot_wider(id_cols = c("protein_id", "peptide_id"), names_from = "sample_id", values_from = "intensity", values_fill = NA)

  if(normalize_intensities) {
    tmp = as.matrix(abundance_table %>% select(-protein_id, -peptide_id))
    tmp = normalize_matrix(tmp, algorithm = "mwmb", group_by_cols = dataset$samples$group[match(colnames(tmp), dataset$samples$sample_id)])
    tmp = normalize_matrix(tmp, algorithm = "modebetween_protein", group_by_rows = abundance_table$protein_id, group_by_cols = dataset$samples$group[match(colnames(tmp), dataset$samples$sample_id)])
    abundance_matrix_nonlog = 2^tmp
    rm(tmp)
  } else {
    abundance_matrix_nonlog = 2^as.matrix(abundance_table %>% select(-protein_id, -peptide_id))
  }


  result = NULL
  for(contr in dataset$contrasts) {
    # set of sample_id for condition 1 and 2
    sid1 = contr$sampleid_condition1
    sid2 = contr$sampleid_condition2
    k1 = length(sid1)
    k2 = length(sid2)

    # determine count threshold
    k1_diff = ifelse(has_kdiff, k_diff, 0L)
    k2_diff = ifelse(has_kdiff, k_diff, 0L)
    if(has_fracdiff) {
      k1_diff = max(k1_diff, as.integer(ceiling(k1 * frac_diff)))
      k2_diff = max(k2_diff, as.integer(ceiling(k2 * frac_diff)))
    }

    # log current contrast and settings, then issue warnings if needed
    append_log(sprintf("contrast: %s - differential detection Nobs_ratio=%.1f, Npep_pass=%d, Nsample1 (max %d) >= %d, Nsample2 (max %d) >= %d", contr$label_contrast, nobs_ratio, npep_pass, k1, k1_diff, k2, k2_diff), type = "info")
    if(k1_diff > k1) {
      append_log(sprintf("not enough samples in condition 1 at given parameters, k1=%d with threshold %d, finding differential detection hits in condition 1 impossible in contrast: %s", k1, k1_diff, contr$label_contrast), type = "warning")
    }
    if(k2_diff > k2) {
      append_log(sprintf("not enough samples in condition 2 at given parameters, k2=%d with threshold %d, finding differential detection hits in condition 2 impossible in contrast: %s", k2, k2_diff, contr$label_contrast), type = "warning")
    }

    # for each condition, count #samples per peptide
    d = count_table %>%
      select(protein_id, peptide_id) %>%
      mutate(
        n1 = matrixStats::rowSums2(count_matrix[ , colnames(count_matrix) %in% sid1, drop = FALSE]),
        n2 = matrixStats::rowSums2(count_matrix[ , colnames(count_matrix) %in% sid2, drop = FALSE])
      ) %>%
      # remove peptides not detected in currently relevant samples (makes the total peptide count per protein specific to this contrast)
      filter(n1 >= 0 | n2 >= 0)

    i = abundance_table %>%
      select(peptide_id) %>%
      mutate(
        int1 = matrixStats::rowSums2(abundance_matrix_nonlog[ , colnames(abundance_matrix_nonlog) %in% sid1, drop = FALSE], na.rm = TRUE),
        int2 = matrixStats::rowSums2(abundance_matrix_nonlog[ , colnames(abundance_matrix_nonlog) %in% sid2, drop = FALSE], na.rm = TRUE)
      )
    i$int1[!is.finite(i$int1)] = 0
    i$int2[!is.finite(i$int2)] = 0

    # classify differentially detected proteins
    x = d %>%
      left_join(i, by = "peptide_id") %>%
      group_by(protein_id) %>%
      summarise(
        npep_total = n(),
        # number of peptides differentially detected in condition 1 as compared to condition 2
        npep_pass1 = sum(n1 - n2 >= k1_diff),
        # vice versa
        npep_pass2 = sum(n2 - n1 >= k2_diff),
        # total number of peptide*sample detection counts
        nobs1 = sum(n1),
        nobs2 = sum(n2),
        # total intensity across peptides*samples
        int1 = sum(int1, na.rm = TRUE),
        int2 = sum(int2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        # % of all possible detections made = number of observed datapoints / (#peptide * #sample)
        fracobs1 = nobs1/(npep_total*k1),
        fracobs2 = nobs2/(npep_total*k2),
        contrast = contr$label,
        `log2fc_nobs2/nobs1` = log2(nobs2/nobs1),
        `log2fc_int2/int1` = log2(int2/int1),
        # test directionally !
        # so more detected in condition 1 and the nobs ratio is in the same direction (and vice versa)
        pass = (npep_pass1 >= npep_pass & nobs1 > nobs_ratio * nobs2 & int1 > int_ratio * int2) |
          (npep_pass2 >= npep_pass & nobs2 > nobs_ratio * nobs1 & int2 > int_ratio * int1)
      ) %>%
      # reorder output
      select(contrast, protein_id, npep_total, npep_pass1, npep_pass2, nobs1, nobs2, fracobs1, fracobs2, `log2fc_nobs2/nobs1`, int1, int2, `log2fc_int2/int1`, pass)

    x$int1[!is.finite(x$int1)] = NA
    x$int2[!is.finite(x$int2)] = NA

    result = bind_rows(result, x)
  }

  return(result)
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
