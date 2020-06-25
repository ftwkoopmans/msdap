
#' Experimental feature: compare the number of peptide detection counts between groups
#'
#' To easily prioritize proteins-of-interest which do not have sufficient abundance values for statistical analysis, but that do have many more detects in one sample group than the other, we provide a simple score based on the count data from the preceding columns.
#'
#' This is an intentionally simplified approach that only serves to rank proteins for further qualitative analysis.
#'
#' First, adjust the peptide detects per sample by the the total number in the sample to account for systematic differences.
#'
#' The zscore metric in diff_detect_zscore_contrast: <contrast>'; A='count_peptides_detected_within_group_A';  B='count_peptides_detected_within_group_A';  result = scale(log2(B+min(Bi)) - log2(A+min(Ai))).    (where scale = subtract the mean, then divide by standard deviation)
#'
#' For all differential detection counts, be careful of over-interpretation and keep differences in sample group size (#replicates) and the total amount of peptides identified in a sample in mind !
#'
#' @param dataset a valid dataset
#'
#' @importFrom data.table dcast
#' @export
differential_detect = function(dataset) {
  # TODO: input validation

  # remove preexisting results
  dataset$dd_proteins = NULL

  column_contrasts = grep("^contrast:", colnames(dataset$samples), ignore.case = T, value = T)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential detection analysis is cancelled", type = "warning")
    return(dataset)
  }


  if(!check_dataset_hascache(dataset)) {
    dataset = cache_filtering_data(dataset)
  }

  # collect count data for output table;
  # per group, count in how many samples a protein was detected
  # per group, count total number of peptide detection events (eg; group with n replicates and m peptides has at most score n*m)
  # note; we use the '_noexclude' fields as we don't take samples into account that were flagged as 'exclude' by the user @ samples input table

  ## 1) collect number of 'detect events' per sample group (outlier samples not taken into account)
  x = dataset$dt_pep_group[ , .(key_protein=key_protein[1],
                                key_group=key_group[1],
                                score = sum(ndetect_noexclude),
                                score_scaled = sum(ndetect_scaled_noexclude)),
                            by=key_protein_group]

  tib_results = tibble()

  ## 2) iterate contrasts
  for (col_contr in column_contrasts) { # col_contr = column_contrasts[2]
    ## 2a) collect score per protein*condition for current contrast
    # unique group * condition pairs for current contrast
    contr_groups = dataset$samples %>%
      select(key_group, condition=!!col_contr) %>%
      filter(condition != 0) %>%
      distinct(key_group, .keep_all = T) %>%
      left_join(dataset$groups, by="key_group")

    condition_size = contr_groups %>% group_by(condition) %>% summarise(size_noexclude = sum(size_noexclude)) %>% arrange(condition)
    stopifnot(condition_size$condition == 1:2)

    # get the sample condition for each key_group. for those groups not in this contrast, result is NA (because of `filter(condition != 0)` above), so remove those. finally, if multi-group; collapse respective score by group
    x_contr = x[, condition := match(key_group, contr_groups$key_group) ][!is.na(condition)][, .(score=sum(score), score_scaled=sum(score_scaled)), by=.(condition, key_protein)]

    ## 2b) differential detect score; score*condition in wide format, then compute z-score
    x_contr_score = data.table::dcast(x_contr, key_protein ~ condition, value.var = "score", fill=0)
    x_contr_score_scaled = data.table::dcast(x_contr, key_protein ~ condition, value.var = "score_scaled", fill=0)
    # assume these wide tables are aligned
    stopifnot(x_contr_score$key_protein == x_contr_score_scaled$key_protein)

    # absolute difference (no scaling)
    zdiff = x_contr_score$`2` - x_contr_score$`1`

    #### log ratio of scaled counts
    ## in the simplified approach, we used the counts as-is (eg; not correcting by the total number of detects in each sample), thus these scores were real count values
    # to handle infinite fold-changes caused by zero's in one condition and not the other, we simply added a + 1
    # zratio_noscale = log2(1 + x_contr_score$`2`) - log2(1 + x_contr_score$`1`)
    ## now in updated approach, we used the count data scaled by sample totals. so to handle the 'no detect' cases on either side, we add the respective minimum scores
    zratio_scaled = log2(x_contr_score_scaled$`2` + min(x_contr_score_scaled$`2`[x_contr_score_scaled$`2` != 0])) -
                    log2(x_contr_score_scaled$`1` + min(x_contr_score_scaled$`1`[x_contr_score_scaled$`1` != 0]))

    # standardize score
    # z_noscale = zratio_noscale - mean(zratio_noscale)
    # z_noscale = z_noscale / sd(z_noscale)
    #
    z_scaled = zratio_scaled - mean(zratio_scaled)
    z_scaled = z_scaled / sd(z_scaled)

    # shortlist candidates; diff in counts is at least 75% of group size  AND  z-score is at least +/- 2
    # z_candidate_noscale = abs(zdiff) >= min(condition_size$size_noexclude)*0.75 & is.finite(z_noscale) & abs(z_noscale) >= 2
    z_candidate_scaled = abs(zdiff) >= min(condition_size$size_noexclude)*0.75 & is.finite(z_scaled) & abs(z_scaled) >= 2

    # append to results
    tib_results = bind_rows(tib_results,
                            tibble(key_protein = x_contr_score$key_protein, diff_detect_zscore=z_scaled, diff_detect_zscore_candidate=z_candidate_scaled, contrast=col_contr)) # , diff_detect_zscore=z_noscale, diff_detect_zscore_candidate=z_candidate_noscale
  }

  # add protein_id and the dd_proteins tibble is done
  tib_results$protein_id = dataset$peptides$protein_id[match(tib_results$key_protein, dataset$peptides$key_protein)]
  # debug; print(tib_results %>% filter(contrast==column_contrasts[2] & diff_detect_zscore_candidate) %>% arrange(desc(abs(diff_detect_zscore))) %>% left_join(dataset$proteins), n=50)
  # debug; print(tib_results %>% filter(contrast==column_contrasts[2] & diff_detect_zscore_corrected_candidate) %>% arrange(desc(abs(diff_detect_zscore_corrected))) %>% left_join(dataset$proteins), n=50)
  dataset$dd_proteins = tib_results
  return(dataset)
}



#' convert the results from differential detection analysis from a long-format tibble to wide-format
#'
#' @param dataset your dataset. if 'dd_proteins' is lacking, result is empty tibble
diffdetect_results_to_wide = function(dataset) {
  if(!is_tibble(dataset$dd_proteins) || nrow(dataset$dd_proteins) == 0) {
    return(tibble())
  }

  ## collect number of 'detect events' per sample group (outlier samples not taken into account)
  x = dataset$dt_pep_group[ , .(key_protein=key_protein[1],
                                key_group=key_group[1],
                                score = sum(ndetect_noexclude)),
                            by=key_protein_group]
  tib_detect_counts = as_tibble(x) %>%
    left_join(dataset$groups %>% select(key_group, group), by="key_group") %>%
    pivot_wider(id_cols = "key_protein", names_from = "group", values_from = "score", values_fill = list(score=0), names_prefix = "count_peptides_detected_within_group_" )

  tib_detect_zscores = dataset$dd_proteins %>%
    pivot_wider(id_cols = "key_protein", names_from = "contrast", values_from = c("diff_detect_zscore", "diff_detect_zscore_candidate")) # , "diff_detect_zscore_corrected", "diff_detect_zscore_corrected_candidate"

  tib_result = full_join(tib_detect_counts, tib_detect_zscores, by="key_protein")
  tib_result$protein_id = dataset$peptides$protein_id[match(tib_result$key_protein, dataset$peptides$key_protein)]
  tib_result$key_protein = NULL

  return(tib_result)
}
