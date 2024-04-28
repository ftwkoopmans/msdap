
#' Summarize DEA and/or differential detection results in a dataset into a table with a single statistic per gene
#'
#' @description
#' For differential detection, the z-scores are converted to p-values by `pnorm(abs(zscore), lower.tail = F)` ,
#' assuming these are normally distributed, and multiple testing correction by FDR is applied.
#' However, this is just an approximation; the differential detection scores are inferior to DEA analyses and should be treated with care.
#' Stringent cutoffs (diffdetect_zscore_threshold of at least 4 ~ 6) are encouraged when using differential detection results.
#'
#' For DEA, the criteria for significance that were previously configured when performing DEA using `dea()` or
#' `analysis_quickstart()` functions are reused here (i.e. cutoffs for adjusted p-value and log2 foldchange).
#'
#' If `return_dea` and `return_diffdetect` are both set to TRUE, the effectsizes from each statistic/result are standardized.
#' This makes both statistics more comparable but the downside is that the effectsizes are no longer the exact effectsizes returned by the DEA model.
#'
#' Note that if your dataset contains results for multiple contrasts and multiple DEA algorithms,
#' these will all be appended into the result table so make sure to filter the results by 'contrast' and
#' 'dea_algorithm' columns where appropriate.
#'
#' @examples \dontrun{
#'   # summarise the DEA results (and ignore the more speculative differential detect data)
#'   x = summarise_stats(dataset, return_dea = TRUE, return_diffdetect = FALSE,
#'                       remove_ambiguous_proteingroups = FALSE)
#'
#'   # first plot all differential detect results as histograms to verify the zscore cutoff
#'   tmp = lapply(plot_differential_detect(dataset), print)
#'   # merge both DEA and differential-detect 'extremes' into a gene-level stats summary
#'   x = summarise_stats(dataset, return_dea = TRUE, return_diffdetect = TRUE,
#'                       diffdetect_zscore_threshold = 5, remove_ambiguous_proteingroups = FALSE)
#' }
#' @param dataset dataset where dea() and/or differential_detect() has been applied
#' @param return_dea boolean, set to TRUE to include DEA results in the stats table that returns 1 value per gene (setting TRUE for both DEA and DD will merge results)
#' @param return_diffdetect analogous to `return_dea`, but setting this to TRUE includes differential detection results
#' @param dea_logfc_as_effectsize optionally, the resulting effectsize column can be populated with standardized foldchange values (effectsize = log2fc / sd(log2fc)).
#' When including differential detection results this'll be a convenient approach to getting 1 standardized distribution of DEA+DD effectsizes that can be used in e.g. GO analyses.
#' While this is unusual, one could e.g. use this for DEA algorithms that apply shrinkage to estimated foldchanges such as MSqRob
#' @param diffdetect_zscore_threshold differential detect z-score cutoff. A typical value would be `diffdetect_zscore_threshold=4` , or 5~6 for stringent filtering.
#' To plot histograms of the respective z-score distributions and inspect potential cutoff values for this relatively arbitrary metric, see below example code
#' @param remove_ambiguous_proteingroups boolean parameter indicating whether proteingroups that contain shared gene symbols (e.g. GRIA1;GRIA2) should be removed
#' @returns \itemize{
#' \item The effectsize column data depends on the user parameters/configuration;
#' When only differential detect is used, this contains the z-scores. When only DEA is used, this contains effectsizes from DEA as-is.
#' However, when DEA and differential detect are combined, the effectsizes from DEA are standardized (divided by std) such that effectsizes from
#' both statistics are integrated into 1 distribution (and thus these can be compared/ranked in downstream analyses).
#' \item 'gene_symbols_or_id' contains the uniprot gene symbols for respective accessions in the proteingroup (protein_id column) as per usual in MS-DAP.
#' \item The 'symbol' column contains the first gene symbol thereof (e.g. in ambiguous groups where gene_symbols_or_id='GRIA1;GRIA2', this'll yield 'GRIA1').
#' \item the 'dea_algorithm' column shows the DEA algorithm that was used. Importantly, even differential detect results that were merged into results for this DEA algorithm have the same dea_algorithm value so you can group/filter by the entire analysis downstream.
#' }
#' @export
summarise_stats = function(dataset, return_dea = TRUE, return_diffdetect = FALSE, dea_logfc_as_effectsize = FALSE, diffdetect_zscore_threshold = 4, remove_ambiguous_proteingroups = FALSE) {
  if(length(return_dea) != 1 || ! return_dea %in% c(TRUE, FALSE)) {
    append_log("return_dea must be single boolean", type = "error")
  }
  if(length(return_diffdetect) != 1 || ! return_diffdetect %in% c(TRUE, FALSE)) {
    append_log("return_diffdetect must be single boolean", type = "error")
  }
  if( ! (return_dea || return_diffdetect) ) {
    append_log("have to set at least one of 'return_dea' and 'return_diffdetect' parameters to TRUE", type = "error")
  }

  all_contrasts = dataset_contrasts(dataset)
  if(length(all_contrasts) == 0) {
    append_log("summarise_stats(); no contrasts in the dataset, returning empty result", type = "warning")
    return(NULL)
  }

  all_algo_dea = NA
  # if the user wants DEA results, find all algorithms used. Halt on unavailable / erronous values
  if(return_dea) {
    if(!("de_proteins" %in% names(dataset) && is.data.frame(dataset$de_proteins) && "dea_algorithm" %in% colnames(dataset$de_proteins))) {
      append_log("DEA results requested by parameter return_dea=TRUE, but no DEA results are available in this dataset. Did you forget to run the analysis_quickstart() or dea() function ?", type = "error")
    }
    all_algo_dea = unique(dataset$de_proteins$dea_algorithm)
    if(length(all_algo_dea) == 0 || anyNA(all_algo_dea) || !is.character(all_algo_dea) || any(all_algo_dea == "")) {
      append_log("DEA results requested by parameter return_dea=TRUE, dataset$de_proteins holds invalid values (missing/NA in 'dea_algorithm' column). Did you forget to run the analysis_quickstart() or dea() function ?", type = "error")
    }
  }


  result = list()
  for(iter_algo_de in all_algo_dea) {
    for(iter_contr in all_contrasts) {
      tmp = summarise_stats__for_contrast(
        dataset,
        return_dea = return_dea,
        return_diffdetect = return_diffdetect,
        contr = iter_contr,
        dea_algorithm = iter_algo_de,
        dea_logfc_as_effectsize = dea_logfc_as_effectsize,
        diffdetect_zscore_threshold = diffdetect_zscore_threshold,
        remove_ambiguous_proteingroups = remove_ambiguous_proteingroups
      )
      if(!is.na(iter_algo_de)) {
        tmp = tmp %>% add_column(dea_algorithm = iter_algo_de, .after = 2)
      }
      result[[length(result) + 1]] = tmp
    }
  }

  bind_rows(result)
}



#' worker function for actual stats summary
#'
#' @param dataset see `summarise_stats()`
#' @param return_dea see `summarise_stats()`
#' @param return_diffdetect see `summarise_stats()`
#' @param contr contrast for which DD and/or DEA was performed (i.e. valid value in dataset$dd_proteins and dataset$de_proteins)
#' @param dea_algorithm DEA algorithm as used in `dea()` upstream
#' @param dea_logfc_as_effectsize see `summarise_stats()`
#' @param diffdetect_zscore_threshold see `summarise_stats()`
#' @param remove_ambiguous_proteingroups see `summarise_stats()`
summarise_stats__for_contrast = function(dataset, return_dea, return_diffdetect, contr, dea_algorithm, dea_logfc_as_effectsize, diffdetect_zscore_threshold, remove_ambiguous_proteingroups) {
  if(!is.list(dataset)) {
    append_log("invalid dataset", type = "error")
  }
  if(!all(c("protein_id", "gene_symbols_or_id") %in% colnames(dataset$proteins)) || !is.character(dataset$proteins$gene_symbols_or_id) || any(dataset$proteins$gene_symbols_or_id == "") ) {
    append_log("requires MS-DAP dataset where the proteins table contains gene symbols (did you forget to import_fasta() ?)", type = "error")
  }
  if(length(return_dea) != 1 || ! return_dea %in% c(TRUE, FALSE)) {
    append_log("return_dea must be single boolean", type = "error")
  }
  if(length(return_diffdetect) != 1 || ! return_diffdetect %in% c(TRUE, FALSE)) {
    append_log("return_diffdetect must be single boolean", type = "error")
  }
  if( ! (return_dea || return_diffdetect) ) {
    append_log("have to set at least one of 'return_dea' and 'return_diffdetect' parameters to TRUE", type = "error")
  }
  if(length(dea_algorithm) != 1 || (!is.na(dea_algorithm) && !is.character(dea_algorithm)) ) {
    append_log("dea_algorithm must be NA (to skip DEA) or a single string", type = "error")
  }
  if(length(contr) != 1 || !is.character(contr)) {
    append_log("contr must be a single string", type = "error")
  }
  if(length(dea_logfc_as_effectsize) != 1 || ! dea_logfc_as_effectsize %in% c(TRUE, FALSE)) {
    append_log("dea_logfc_as_effectsize must be single boolean", type = "error")
  }
  if(length(diffdetect_zscore_threshold) != 1 || !is.numeric(diffdetect_zscore_threshold) || !is.finite(diffdetect_zscore_threshold) || diffdetect_zscore_threshold < 0) {
    append_log("diffdetect_zscore_threshold must be a single positive numeric value", type = "error")
  }
  if(length(remove_ambiguous_proteingroups) != 1 || ! remove_ambiguous_proteingroups %in% c(TRUE, FALSE)) {
    append_log("remove_ambiguous_proteingroups must be single boolean", type = "error")
  }

  param_dea_algorithm = dea_algorithm

  log = dd_type = ""
  result = NULL
  dea = dd = data.frame()
  use_dea = return_dea && !is.na(param_dea_algorithm)
  use_dd = return_diffdetect
  if(!use_dea && !use_dd) {
    append_log("invalid paramters, both DEA and DD are disabled; dea_algorithm=NA and diffdetect_zscore_threshold=NA", type = "error")
  }


  # user wants to use DEA data, double-check dataset contains DEA results for selected contrast and dea_algorithm
  if(use_dd) {
    if(!"dd_proteins" %in% names(dataset) || !all(c("protein_id", "contrast", "zscore") %in% colnames(dataset$dd_proteins)) ) {
      append_log("if 'return_diffdetect' parameter is set to TRUE, an MS-DAP dataset with differential detection z-scores is required that were computed using the `differential_detect()` function with parameter 'return_wide_format = FALSE'", type = "error")
    }
    if( ! contr %in% dataset$dd_proteins$contrast) {
      append_log(sprintf("selected contrast '%s' is not available in the differential detection results; please re-run differential_detect() to ensure data is available for all contrasts specified for this dataset", contr), type = "error")
    }

    dd = dataset$dd_proteins %>% filter(contrast == contr & is.finite(zscore))

    if(nrow(dd) > 0) {
      # results might contain both detect and quant, or only one of the other
      dd_type = unique(dd$type)
      if("detect" %in% dd_type) {
        dd_type = "detect"
      } else {
        dd_type = dd_type[1]
      }

      # subset the DD results for the selected contrast @ type and z-score cutoff
      dd = dd %>%
        filter(type == dd_type) %>%
        mutate(
          zscore_pvalue = stats::pnorm(abs(zscore), lower.tail = F),
          zscore_pvalue_adjust = p.adjust(zscore_pvalue, method = "BH"),
          signif = abs(zscore) >= diffdetect_zscore_threshold
        ) %>%
        arrange(desc(abs(zscore))) %>%
        select(protein_id, peptides_used_for_dd = npep_max, log2fc_dd = log2fc, effectsize_dd = zscore, pvalue_dd = zscore_pvalue, pvalue_adjust_dd = zscore_pvalue_adjust, signif_dd = signif)
    }

    if(nrow(dd) == 0) {
      log = sprintf("%sconfigured to use differential detect, but no results match user parameters. ", log)
    }
  }



  # user wants to use DEA data, double-check dataset contains DEA results for selected contrast and dea_algorithm
  if(use_dea) {
    if(!"de_proteins" %in% names(dataset) || !all(c("protein_id", "contrast", "dea_algorithm", "pvalue", "qvalue", "foldchange.log2", "effectsize", "peptides_used_for_dea", "signif") %in% colnames(dataset$de_proteins)) ) {
      append_log("DEA results are requested but the provided MS-DAP dataset has none. Did you forget to run the analysis_quickstart() or dea() function ?", type = "error")
    }
    if( ! contr %in% dataset$de_proteins$contrast) {
      append_log(sprintf("selected contrast '%s' is not available in the DEA results; please re-run dea() to ensure data is available for all contrasts specified for this dataset", contr), type = "error")
    }

    dea = dataset$de_proteins %>%
      # importantly, filter for contrast, DEA algorithm and #peptides
      filter(dea_algorithm == param_dea_algorithm & is.finite(qvalue) & contrast == contr) %>%
      mutate(effectsize_dea = effectsize,
             effectsize_dea_std = effectsize_dea / sd(effectsize_dea, na.rm=T),
             log2fc_dea = foldchange.log2,
             log2fc_dea_std = log2fc_dea / sd(log2fc_dea, na.rm=T),
             pvalue_dea = pvalue,
             pvalue_adjust_dea = qvalue
      ) %>%
      select(protein_id, peptides_used_for_dea, log2fc_dea, log2fc_dea_std, effectsize_dea, effectsize_dea_std, pvalue_dea, pvalue_adjust_dea, signif_dea = signif) %>%
      arrange(pvalue_dea)

    if(dea_logfc_as_effectsize) {
      dea$effectsize_dea = dea$log2fc_dea
      dea$effectsize_dea_std = dea$log2fc_dea_std
    }

    if(nrow(dea) == 0) {
      log = sprintf("%sconfigured to use DEA, but no results match user parameters (e.g. no data for combination of 'dea_algorithm' and 'contrast'). ", log)
    }
  }


  if(nrow(dd) == 0 && nrow(dea) == 0) {
    append_log("filtering settings did not yield any proteins", type = "info")
    return(NULL)
  }


  # only DEA
  if(nrow(dea) > 0 && (nrow(dd) == 0 || !any(dd$signif_dd))) {
    result = dea %>% mutate(score_type = "dea") %>% select(protein_id, score_type, peptide_count = peptides_used_for_dea, log2fc = log2fc_dea, effectsize = effectsize_dea, pvalue = pvalue_dea, pvalue_adjust = pvalue_adjust_dea, signif = signif_dea)
    log = sprintf("%sreturning only DEA results from %s, %d proteins.", log, param_dea_algorithm, nrow(result))
  }

  # only DD
  if(nrow(dea) == 0 && nrow(dd) > 0) {
    result = dd %>% mutate(score_type = "dd") %>% select(protein_id, score_type, peptide_count = peptides_used_for_dd, log2fc = log2fc_dd, effectsize = effectsize_dd, pvalue = pvalue_dd, pvalue_adjust = pvalue_adjust_dd, signif = signif_dd)
    log = sprintf("%sreturning only differential detect results, %d proteins.", log, nrow(result))
  }

  # merge both, with DEA as leading and only hits from DD that make it past the zscore threshold (function parameter)
  if(nrow(dea) > 0 && nrow(dd) > 0 && any(dd$signif_dd)) {

    result = full_join(dea, dd %>% filter(signif_dd), by = "protein_id") %>%
      # default result values; use DEA. note that we use the standardized effect-size so we can merge with DD
      mutate(peptide_count = peptides_used_for_dea, log2fc = log2fc_dea, effectsize = effectsize_dea_std, pvalue = pvalue_dea, pvalue_adjust = pvalue_adjust_dea,
             score_type = ifelse(is.na(pvalue_dea), NA, "dea"))

    # overwrite DEA statistics where Differential Detect (DD) has stronger result
    # (but note that we only consider the set of 'strong' DD scores, e.g. we don't overwrite poor DEA p-values with slightly better DD results)
    rows_na        = !is.finite(result$pvalue_dea) & is.finite(result$pvalue_dd)
    rows_stronger  = is.finite(result$pvalue_dea) & is.finite(result$pvalue_dd) & result$pvalue_dea > result$pvalue_dd
    rows = rows_na | rows_stronger
    result$score_type[rows]    = ifelse(is.na(result$score_type[rows]), "dd", "dea,dd")
    result$log2fc[rows]        = result$log2fc_dd[rows]
    result$effectsize[rows]    = result$effectsize_dd[rows]
    result$pvalue[rows]        = result$pvalue_dd[rows]
    result$pvalue_adjust[rows] = result$pvalue_adjust_dd[rows]
    result$peptide_count[rows] = result$peptides_used_for_dd[rows]
    # significant; any protein that remains after differential-detect cutoff or signif @ DEA
    result$signif = is.finite(result$effectsize_dd) | result$signif_dea %in% TRUE # `%in% TRUE` also deals with NA values

    result = result %>%
      filter(is.finite(pvalue)) %>%
      select(protein_id, score_type, peptide_count, log2fc, effectsize_dea, effectsize, pvalue, pvalue_adjust, signif)

    log = sprintf("%smerging DEA results from %s and differential detect (score type '%s'), the latter contributed to %d / %d proteingroups (%d without DEA, %d stronger than DEA).",
                  log, param_dea_algorithm, dd_type, sum(rows_na | rows_stronger), length(rows_na), sum(rows_na), sum(rows_stronger))
  }


  if(!is.data.frame(result) || nrow(result) == 0) {
    append_log("filtering settings did not yield any proteins", type = "info")
    return(NULL)
  }

  result = result %>%
    # add (leading) gene symbols
    left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by = "protein_id") %>%
    mutate(symbol = gsub(";.*", "", gene_symbols_or_id)) %>%
    # add contrast
    add_column(contrast = contr, .after = 1) %>%
    arrange(desc(signif), pvalue)


  # finally, retain 1 row per gene
  # this crucially depends on upstream sorting of the table
  if(remove_ambiguous_proteingroups) {
    # remove all proteingroups with multiple genes, then retain unique genes
    result = result %>% filter(!grepl(";", gene_symbols_or_id, fixed = TRUE)) %>% distinct(gene_symbols_or_id, .keep_all = T)
  } else {
    # retain all unique shared gene symbols
    result = result %>% distinct(gene_symbols_or_id, .keep_all = T)
  }

  append_log(sprintf(
    "%s;  %s  %d/%d(%.1f%%) gene symbols in the final result table are foreground.",
    contr, log,
    sum(result$signif), nrow(result), sum(result$signif)/nrow(result) * 100
  ), type = "info")

  return(result)
}
