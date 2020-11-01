
#' Compare the number of peptide detection counts between groups
#'
#' Compute a z-score for each protein based on the total number of detected peptides per sample group.
#'
#' To easily prioritize proteins-of-interest, some which may not have sufficient abundance values for differential expression analysis but that do have many more detects in one sample group than the other, we provide a simple score based on identification count data.
#'
#' This is a simplified approach that is intended to rank proteins for further qualitative analysis, be careful of over-interpretation and keep differences in sample group size (#replicates) and the absolute amount of peptides identified in a sample in mind !
#'
#' @param dataset a valid dataset
#' @param min_samples_observed minimum number of samples where a protein should be observed at least once by any of its peptides (in either group) when comparing a contrast of group A vs B
#'
#' @importFrom data.table dcast
#' @importFrom tidyr pivot_wider everything
#' @importFrom matrixStats rowSums2
#' @export
differential_detect = function(dataset, min_samples_observed = 3) { #, threshold_fraction_detect = 0.5
  stopifnot(length(min_samples_observed) == 1 && !is.na(min_samples_observed) && is.numeric(min_samples_observed) && min_samples_observed > 0)

  if("fraction" %in% colnames(dataset$samples)) {
    append_log("you have to merge fractionated samples first. Typically, you want to simply use the analysis_quickstart() to process the dataset first. Advanced users may call merge_fractionated_samples(dataset) manually.", type = "error")
  }

  # remove preexisting results
  dataset$dd_proteins = NULL

  column_contrasts = dataset_contrasts(dataset)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential detection analysis is cancelled", type = "warning")
    return(dataset)
  }


  if(!check_dataset_hascache(dataset)) {
    dataset = cache_filtering_data(dataset)
  }

  append_log(sprintf("differential detection analysis: min_samples_observed=%d", min_samples_observed), type = "info")


  # collect count data for output table;
  # per group, count in how many samples a protein was detected
  # per group, count total number of peptide detection events (eg; group with n replicates and m peptides has at most score n*m)
  # note; we use the '_noexclude' fields as we don't take samples into account that were flagged as 'exclude' by the user @ samples input table

  ## 1A) collect number of 'detect events' per sample group (outlier samples not taken into account)
  count_data = list(
    detect = dataset$dt_pep_group[ , .(key_protein=key_protein[1],
                                       key_group=key_group[1],
                                       score = sum(ndetect_noexclude),
                                       score_scaled = sum(ndetect_scaled_noexclude)),
                                   by=key_protein_group]
  )

  if(!is_dia_dataset(dataset)) {
    # analogous, but for quantification counts (have some abundance / intensity, may or may not be 'detected')
    count_data[["quant"]] = quant = dataset$dt_pep_group[ , .(key_protein=key_protein[1],
                                                              key_group=key_group[1],
                                                              score = sum(nquant_noexclude),
                                                              score_scaled = sum(nquant_scaled_noexclude)),
                                                          by=key_protein_group]
  }


  ## 1B) collect number of samples where each protein was detected at least once
  ## this is useful downstream because we might not be interested in a protein that has 5 peptides where all peptides were detected in 1 out of 5 WT samples and never in knockout (absolute detect count difference is 5, but all of it relies on 1 sample)
  # in below tables, all samples flagged 'exclude' are disregarded
  key_sample__not_excluded = dataset$samples %>% filter(exclude != TRUE) %>% pull(key_sample)
  # note; proteins that have some peptide intensity but are not detected (eg; match-between-runs in DDA) are recognised as 'detect==FALSE' in this first summary table, proteins that are totally absent have no entry/row in this sample*protein table
  tib_count_protein_by_sample = dataset$peptides %>% filter(key_sample %in% key_sample__not_excluded) %>% group_by(key_sample, key_protein) %>% summarise(detect=any(detect)) %>% ungroup()
  # gather key_group from samples table (instead of peptides because that table is huge which makes lookups slower)
  tib_count_protein_by_sample$key_group = dataset$samples$key_group[match(tib_count_protein_by_sample$key_sample, dataset$samples$key_sample)]
  # any peptide of protein p was found in sample s -->> count number of samples where p was detected per group
  tib_count_protein_by_group = tib_count_protein_by_sample %>% group_by(key_group, key_protein) %>% summarise(ndetect = sum(detect), nquant = n()) %>% ungroup()

  # cast to wide table for downstream efficiency (so we can use matrix operations instead of long-format table group_by and summary 2*ncontrast times)
  tibw = tib_count_protein_by_group %>% pivot_wider(id_cols = "key_protein", names_from = "key_group", values_from = "ndetect")
  mat_protein_by_group__key_protein = tibw$key_protein
  mat_protein_by_group_detect = as_matrix_except_first_column(tibw)
  mat_protein_by_group_detect[!is.finite(mat_protein_by_group_detect)] = 0
  rm(tibw)

  if(!is_dia_dataset(dataset)) {
    mat_protein_by_group_quant = as_matrix_except_first_column(tib_count_protein_by_group %>% pivot_wider(id_cols = "key_protein", names_from = "key_group", values_from = "nquant"))
    mat_protein_by_group_quant[!is.finite(mat_protein_by_group_quant)] = 0
    #
    stopifnot(rownames(mat_protein_by_group_detect) == rownames(mat_protein_by_group_quant))
  }


  tib_results = tibble()
  for(type_detect_quant in names(count_data)) { # type_detect_quant=names(count_data)[1]
    x = count_data[[type_detect_quant]]

    ## 2) iterate contrasts
    for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
      ## 2a) collect score per protein*condition for current contrast
      # unique group * condition pairs for current contrast
      contr_groups = dataset$samples %>%
        select(key_group, condition=!!col_contr) %>%
        filter(condition != 0) %>%
        distinct(key_group, .keep_all = T) %>%
        left_join(dataset$groups, by="key_group") %>%
        arrange(condition)

      condition_size = contr_groups %>% group_by(condition) %>% summarise(size_noexclude = sum(size_noexclude)) %>% arrange(condition)
      stopifnot(condition_size$condition == 1:2)

      # get the sample condition for each key_group. for those groups not in this contrast, result is NA (because of `filter(condition != 0)` above), so remove those. finally, if multi-group; collapse respective score by group
      x_contr = x[, condition := contr_groups$condition[match(key_group, contr_groups$key_group)] ][!is.na(condition)][, .(score=sum(score), score_scaled=sum(score_scaled)), by=.(condition, key_protein)]


      ## count the number of samples where each protein was observed per condition
      m = mat_protein_by_group_detect
      if(type_detect_quant == "quant") {
        m = mat_protein_by_group_quant
      }
      # !! as.character is pivotal, R matrices with numeric row and col names are a nightmare
      g1 = contr_groups %>% filter(condition == 1) %>% pull(key_group) %>% as.character
      g2 = contr_groups %>% filter(condition == 2) %>% pull(key_group) %>% as.character
      # for each condition, the total count of 'samples where protein p was observed' is the sum of respective sample groups
      # next, we take the maximum over both conditions for each protein. This is matched against user filtering criteria (eg; proteins only seen in 2 out of 5 samples _at most_ are not interesting regardless of z-score)
      nsample_observed_max = tibble(key_protein = mat_protein_by_group__key_protein, # these are the rownames in m, but as integers (cached earlier)
                                    nmax = pmax(matrixStats::rowSums2(m[,as.character(g1),drop=F]), matrixStats::rowSums2(m[,as.character(g2),drop=F])) ) # redundantly enforce as.character  (not needed, but safety first around such bug-prone code)

      ## reference code, in case we don't have wide-format tibble (this is much slower)
      # # - by inner-joining the groups used in current contrast, we subset the total count table for the sample groups used in current contrast
      # # - then summarise counts at condition level by summation
      # # - finally, we can classify which proteins pass the user-defined selection criterium
      # tib_nsamples = tib_count_protein_by_group %>% inner_join(contr_groups %>% select(key_group, condition), by = "key_group") %>%
      #   group_by(key_protein, condition) %>% summarise(ndetect = sum(ndetect), nquant = sum(nquant)) %>%
      #   group_by(key_protein) %>% summarise(ndetect = max(ndetect), nquant = max(nquant)) %>%
      #   ungroup()
      # if(type_detect_quant == "detect") {
      #   tib_nsamples$candidate = tib_nsamples$ndetect >= min_samples_observed
      # } else {
      #   tib_nsamples$candidate = tib_nsamples$nquant >= min_samples_observed
      # }


      x_contr_score_scaled = as_tibble(data.table::dcast(x_contr, key_protein ~ condition, value.var = "score_scaled", fill=0)) %>%
        # add 'max replicate samples over conditions' to the count tibble used for z-scores
        left_join(nsample_observed_max, by = "key_protein")

      # only proteins that pass user-specified filtering criteria; at least detected in N samples in either condition 1 or 2
      x_contr_score_scaled = x_contr_score_scaled %>% filter(is.finite(nmax) & nmax >= min_samples_observed)

      # log2 scaled foldchange, padding each condition with minimum non-zero score (to deal with zero cases)
      zratio_scaled = log2(x_contr_score_scaled$`2` + min(x_contr_score_scaled$`2`[x_contr_score_scaled$`2` != 0])) -
        log2(x_contr_score_scaled$`1` + min(x_contr_score_scaled$`1`[x_contr_score_scaled$`1` != 0]))

      # standardize score
      z_scaled = zratio_scaled - mean(zratio_scaled, na.rm = T)
      z_scaled = z_scaled / sd(z_scaled, na.rm = T)

      tib_results = bind_rows(tib_results,
                              tibble(key_protein = x_contr_score_scaled$key_protein, zscore_count=z_scaled, contrast=col_contr, type=type_detect_quant)) # , zscore_count=z_noscale, zscore_count_candidate=z_candidate_noscale

      ################ v1
      # ## differential detect score; score*condition in wide format, then compute z-score
      # x_contr_score = data.table::dcast(x_contr, key_protein ~ condition, value.var = "score", fill=0)
      # x_contr_score_scaled = data.table::dcast(x_contr, key_protein ~ condition, value.var = "score_scaled", fill=0)
      # # enforce table sorting such that both wide-format tables are aligned at all times
      # x_contr_score = x_contr_score[ , c("key_protein", "1", "2")]
      # x_contr_score_scaled = x_contr_score_scaled[match(x_contr_score_scaled$key_protein, x_contr_score$key_protein) , c("key_protein", "1", "2")] # protein order exact same as previous (should be by default but enforce anyway)
      # stopifnot(!is.na(x_contr_score_scaled$key_protein))
      # # QC plot to check scaled variant is strongly correlated;
      # # plot(x_contr_score$`1`, x_contr_score_scaled$`1`)
      # # plot(x_contr_score$`2`, x_contr_score_scaled$`2`)
      #
      # ## add 'max replicate samples over conditions' to the count tibble used for z-scores
      # # index in the `nsample_observed_max` table for each row in `x_contr_score`
      # i = match(x_contr_score$key_protein, nsample_observed_max$key_protein)
      # i_isna = is.na(i)
      # x_contr_score_scaled$nmax = nsample_observed_max$nmax[i] # x_contr_score$nmax =
      # x_contr_score$nmax[i_isna] = 0
      # x_contr_score_scaled$nmax[i_isna] = 0
      #
      #
      # # # absolute difference (no scaling)
      # # zdiff = x_contr_score$`2` - x_contr_score$`1`
      #
      # #### log ratio of scaled counts
      # ## in the simplified approach, we used the counts as-is (eg; not correcting by the total number of detects in each sample), thus these scores were real count values
      # # to handle infinite fold-changes caused by zero's in one condition and not the other, we simply added a + 1
      # # zratio_noscale = log2(1 + x_contr_score$`2`) - log2(1 + x_contr_score$`1`)
      # ## now in updated approach, we used the count data scaled by sample totals. so to handle the 'no detect' cases on either side, we add the respective minimum scores
      #
      #
      # zratio_scaled = log2(x_contr_score_scaled$`2` + min(x_contr_score_scaled$`2`[x_contr_score_scaled$`2` != 0])) -
      #                 log2(x_contr_score_scaled$`1` + min(x_contr_score_scaled$`1`[x_contr_score_scaled$`1` != 0]))
      #
      # # standardize score
      # # z_noscale = zratio_noscale - mean(zratio_noscale)
      # # z_noscale = z_noscale / sd(z_noscale)
      # #
      # z_scaled = zratio_scaled - mean(zratio_scaled, na.rm = T)
      # z_scaled = z_scaled / sd(z_scaled, na.rm = T)
      #
      # # shortlist candidates; diff in counts is at least 75% of group size  AND  z-score is at least +/- 2
      # z_candidate_scaled = x_contr_score_scaled$nmax >= min_samples_observed  &  is.finite(z_scaled) & abs(z_scaled) >= 2
      # # z_candidate_noscale = abs(zdiff) >= min(condition_size$size_noexclude)*0.75 & is.finite(z_noscale) & abs(z_noscale) >= 2
      # # z_candidate_scaled = abs(zdiff) >= min(condition_size$size_noexclude)*threshold_fraction_detect & is.finite(z_scaled) & abs(z_scaled) >= 2
      #
      # # append to results
      # tib_results = bind_rows(tib_results,
      #                         tibble(key_protein = x_contr_score$key_protein, zscore_count=z_scaled, zscore_count_candidate=z_candidate_scaled, contrast=col_contr, type=type_detect_quant)) # , zscore_count=z_noscale, zscore_count_candidate=z_candidate_noscale
      ################ v1
    }
  }

  tib = tidyr::pivot_wider(tib_results, id_cols = c("key_protein", "contrast"), names_from = "type", names_prefix = "zscore_count_", values_from = "zscore_count") # , "zscore_count_candidate"

  # finally, replace the "key" we use internally (faster to match integers all the time) with the actual protein_id
  tib$protein_id = dataset$peptides$protein_id[match(tib$key_protein, dataset$peptides$key_protein)]

  # reorder column names
  dataset$dd_proteins = tib %>% select(protein_id, key_protein, contrast, dplyr::starts_with("zscore"))
  return(dataset)

  # debug; print(tib %>% filter(contrast==column_contrasts[2] & zscore_count_candidate) %>% arrange(desc(abs(zscore_count))) %>% left_join(dataset$proteins), n=50)
  # debug; print(tib %>% filter(contrast==column_contrasts[2] & zscore_count_corrected_candidate) %>% arrange(desc(abs(zscore_count_corrected))) %>% left_join(dataset$proteins), n=50)
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
  x = dataset$dt_pep_group[ , .(key_protein=key_protein[1],
                                key_group=key_group[1],
                                count_peptides_detected_within_group = sum(ndetect_noexclude),
                                count_peptides_quantified_within_group = sum(nquant_noexclude)),
                            by=key_protein_group]

  tib_counts = as_tibble(x) %>%
    left_join(dataset$groups %>% select(key_group, group), by="key_group") %>%
    tidyr::pivot_wider(id_cols = "key_protein", names_from = "group", names_sep = "__", values_from = c("count_peptides_detected_within_group", "count_peptides_quantified_within_group"), values_fill = list(count_peptides_detected_within_group=0, count_peptides_quantified_within_group=0))

  tib_zscores = dataset$dd_proteins %>%
    tidyr::pivot_wider(id_cols = "key_protein", names_from = "contrast", values_from = intersect(c("zscore_count_detect", "zscore_count_quant"), colnames(dataset$dd_proteins)) )


  tib_result = full_join(tib_counts, tib_zscores, by="key_protein")
  tib_result$protein_id = dataset$peptides$protein_id[match(tib_result$key_protein, dataset$peptides$key_protein)]
  tib_result$key_protein = NULL

  return(tib_result %>% select(protein_id, tidyr::everything()))
}
