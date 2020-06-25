
#' returns named variable variable (label=column_name) indicating which type of intensity data is used
#' @param peptides todo
#' @param contr_lbl todo
get_column_intensity = function(peptides, contr_lbl = NA) {
  ref = c("filter by contrast" = paste0("intensity_", contr_lbl),
          "global data filter" = "intensity_all_group",
          "custom filtered/normalized" = "intensity_norm",
          "input data as-is" = "intensity")
  return(ref[ref %in% colnames(peptides)][1])
}



#' Differential expression analysis
#'
#' @param dataset your dataset
#' @param qval_signif threshold for significance of adjusted p-values. default: 0.05
#' @param fc_signif threshold for significance of log2 foldchanges. Set to zero to disregard, a positive value to apply a cutoff to absolute foldchanges or use bootstrap analyses to infer a suitable foldchange threshold by providing either NA or a negative value. default: 0
#' @param algo_de algorithms for differential expression analysis. options: ebayes, msempire, msqrob (to run multiple, provide an array)
#' @param norm_modebetween_protein_eset normalise protein-level data matrix used for DEA by 'modebetween' to correct for 'imbalance' between conditions introduced in peptide-to-protein rollup. Note that some algo_de, such as msqrob or msempire, directly operate on the peptide-level data thus are unaffected by this setting. Only useful for statistical models that first roll-up to protein-level, such as eBayes. default:FALSE
#' @param output_dir_for_eset optionally, provide an output directory where the expressionset objects should be stored. Only useful if you're doing downstream data analysis that required this data
#' @export
dea = function(dataset, qval_signif = 0.05, fc_signif=0, algo_de = c("ebayes"), norm_modebetween_protein_eset = FALSE, output_dir_for_eset = "") {
  ### input validation
  if (length(qval_signif) != 1 || !is.finite(qval_signif) || qval_signif <= 0) {
    append_log("q-value threshold must be a single numerical value above 0 (parameter qval_signif)", type = "error")
  }

  if (length(fc_signif) != 1 || (!is.na(fc_signif) && !is.finite(fc_signif))) {
    append_log("log2 foldchange threshold must be a single numerical value or NA (parameter fc_signif)", type = "error")
  }

  # trigger automatic estimation of FC by either NA or a negative value (to convenience the user, no functional implications)
  if(!is.na(fc_signif) && fc_signif < 0) {
    fc_signif = NA
  }

  algo_de = unique(algo_de)
  if (length(algo_de) == 0) {
    append_log("no algorithms have been defined (parameter algo_de), differential expression analysis is cancelled", type = "warning")
    return(dataset)
  }

  if (length(norm_modebetween_protein_eset) != 1 || is.na(norm_modebetween_protein_eset) || !is.logical(norm_modebetween_protein_eset)) {
    append_log("parameter norm_modebetween_protein_eset must be a boolean value (TRUE or FALSE)", type = "error")
  }

  # valid DEA options are those hardcoded, or pre-existing functions
  global_func = ls(envir=.GlobalEnv)
  algo_de_invalid = setdiff(algo_de, c("ebayes", "msempire", "msqrob", "msqrobsum", global_func))
  if (length(algo_de_invalid) > 0) {
    append_log(paste("invalid options for algo_de:", paste(algo_de_invalid, collapse=", ")), type = "error")
  }

  column_contrasts = grep("^contrast:", colnames(dataset$samples), ignore.case = T, value = T)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential expression analysis is cancelled", type = "warning")
    return(dataset)
  }


  ### data checks out
  # remove preexisting results
  dataset$de_proteins = NULL
  # init result tibble
  result_stats = tibble()

  ### iterate contrasts
  for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
    append_log(paste("differential abundance analysis for", col_contr), type = "info")

    # returns named variable variable (label=column_name) indicating which type of intensity data is used
    col_contr_intensity = get_column_intensity(dataset$peptides, col_contr)
    append_log(paste("using data from peptide filter:", names(col_contr_intensity)), type = "info")

    # if there are no intensity values, there is nothing to do
    if(sum(!is.na(dataset$peptides %>% pull(!!col_contr_intensity))) < 2) {
      append_log("no peptides with an intensity value are available (perhaps filtering was too stringent?)", type = "warning")
      next
    }

    # select the current contrast as column 'condition' and remove irrelevant samples (when defining our contrasts in upstream code, irrelevant samples were designated a 0)
    # ! sorting by 'condition' here is important, as it aligns the samples by the groups specified by the user. if we don't enforce this, an intended WT-vs-KO might actually be analyzed as KO-vs-WT depending on the input data's sample ordering
    samples_for_contrast = dataset$samples %>%
      select(sample_id, shortname, group, condition = !!col_contr) %>%
      filter(condition != 0) %>%
      arrange(condition)

    if(!all(samples_for_contrast$sample_id %in% dataset$peptides$sample_id)) {
      append_log("all samples from this contrast to be in the peptides tibble", type = "error")
    }

    ## convert our long-format peptide table to a peptide- and protein-level ExpressionSet
    peptides_for_contrast = dataset$peptides %>%
      select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, intensity=!!as.character(col_contr_intensity)) %>%
      filter(sample_id %in% samples_for_contrast$sample_id & is.finite(intensity))

    eset_peptides = tibble_as_eset(peptides_for_contrast, dataset$proteins, samples_for_contrast)
    eset_proteins = eset_from_peptides_to_proteins(eset_peptides)

    if(norm_modebetween_protein_eset) {
      append_log("normalizing protein-level data for DEA by 'modebetween' (doesn't affect peptide-level models)", type = "info")
      Biobase::exprs(eset_proteins) = normalize_matrix(Biobase::exprs(eset_proteins), algorithm = "modebetween", mask_sample_groups = Biobase::pData(eset_proteins)$condition)
    }

    # if a directory for file storage was provided, store eset in a .RData file
    if(length(output_dir_for_eset) == 1 && !is.na(output_dir_for_eset) && nchar(output_dir_for_eset)>0 && dir.exists(output_dir_for_eset)) {
      save(eset_peptides, file=paste0(output_dir_for_eset, "/ExpressionSet_peptides_", gsub("\\W+", " ", col_contr), ".RData"), compress = T)
      save(eset_proteins, file=paste0(output_dir_for_eset, "/ExpressionSet_proteins_", gsub("\\W+", " ", col_contr), ".RData"), compress = T)
    }

    contr_fc_signif = fc_signif
    if(is.na(fc_signif)) {
      contr_fc_signif = dea_protein_background_foldchange_limits(eset_proteins)
      append_log(sprintf("log2 foldchange threshold estimated by bootstrap analysis: %.3f", contr_fc_signif), type = "info")
    }

    # DE statistics for all requested algorithms
    tib = tibble()
    for(alg in algo_de) {
      alg_matched = F
      if (alg == "ebayes") {
        # DEBUG <<- eset_proteins
        tib = bind_rows(tib, de_interface_ebayes(eset_proteins=eset_proteins, input_intensities_are_log2 = T))
        alg_matched = T
      }
      if (alg == "msempire") {
        tib = bind_rows(tib, de_msempire(eset_peptides, input_intensities_are_log2 = T) %>% add_column(algo_de = "msempire"))
        alg_matched = T
      }
      if (alg == "msqrob") {
        tib = bind_rows(tib, de_msqrobsum_msqrob(eset_peptides, use_peptide_model = T, input_intensities_are_log2 = T) %>% add_column(algo_de = "msqrob"))
        alg_matched = T
      }
      if (alg == "msqrobsum") {
        tib = bind_rows(tib, de_msqrobsum_msqrob(eset_peptides, use_peptide_model = F, input_intensities_are_log2 = T) %>% add_column(algo_de = "msqrobsum"))
        alg_matched = T
      }
      # for non-hardcoded functions, we call the function requested as a user parameter and pass all available data
      if(!alg_matched) {
        alg_fun = match.fun(alg)
        alg_result = alg_fun(peptides=peptides_for_contrast, samples=samples_for_contrast, eset_peptides=eset_peptides, eset_proteins=eset_proteins, input_intensities_are_log2 = T)
        # validation checks on expected output, to facilitate debugging/feedback for custom implementations
        if(!is_tibble(alg_result) || nrow(alg_result) == 0 || !all(c("protein_id", "pvalue", "qvalue", "foldchange.log2", "algo_de") %in% colnames(alg_result))) {
          append_log(sprintf("provided custom function for differential abundance analysis '%s' must return a non-empty tibble with the columns protein_id|pvalue|qvalue|foldchange.log2|algo_de", alg), type = "error")
        }
        alg_result_name = unique(alg_result$algo_de)
        if(length(alg_result_name) != 1 || !is.character(alg_result_name) || nchar(alg_result_name) < 2 || alg_result_name %in% c("ebayes", "msempire", "msqrob", "msqrobsum", "combined")) {
          append_log(sprintf("provided custom function for differential abundance analysis '%s' contains invalid values in algo_de column (must be a single non-empty character string uniquely indicating the name/label of your method. cannot be either of ebayes|msempire|msqrob|msqrobsum|combined)", alg), type = "error")
        }
        if(!all(!is.na(alg_result$protein_id) & alg_result$protein_id %in% peptides_for_contrast$protein_id)) {
          append_log(sprintf("provided custom function for differential abundance analysis '%s' contains invalid values in protein_id column (either NA or not found in provided data structures)", alg), type = "error")
        }
        if(!all( (is.na(alg_result$pvalue) | is.numeric(alg_result$pvalue)) &
                 (is.na(alg_result$qvalue) | is.numeric(alg_result$qvalue)) &
                 (is.na(alg_result$foldchange.log2) | is.numeric(alg_result$foldchange.log2)) ))  {
          append_log(sprintf("provided custom function for differential abundance analysis '%s' contains invalid values in pvalue|qvalue|foldchange.log2 columns (can only be NA or numeric)", alg), type = "error")
        }
        # finally, concatenate results
        tib = bind_rows(tib, alg_result)
      }
    }

    if(nrow(tib) > 0) {
      # count peptides per protein
      prot_pep_count = peptides_for_contrast %>% distinct(protein_id, peptide_id) %>% count(protein_id, name = "peptides_used_for_dea")
      # add a column with significance flag
      s = is.finite(tib$qvalue) & tib$qvalue <= qval_signif
      if(is.finite(contr_fc_signif) & contr_fc_signif > 0) {
        s = s & is.finite(tib$foldchange.log2) & abs(tib$foldchange.log2) >= contr_fc_signif
      }
      result_stats = bind_rows(result_stats, tib %>% left_join(prot_pep_count, by="protein_id") %>% add_column(signif = s, signif_threshold_qvalue = qval_signif, signif_threshold_log2fc = contr_fc_signif, contrast = col_contr, .after = "qvalue"))
    }
  }

  dataset$de_proteins = result_stats
  return(dataset)
}



#' estimate a threshold for 'significant' foldchanges from N permutations of sample-to-condition assignments
#'
#' Permutations of sample labels within a group are disregarded as these have no effect on the between-group foldchange, only unique combinations of swapping samples between conditions A and B are considered
#'
#' This is somewhat similar to the method described by Hafemeister and Satija at https://doi.org/10.1186/s13059-019-1874-1
#'
#' M&M quote: "A random background distribution of mean differences was generated by randomly choosing 1000 genes and permuting the group labels. Significance thresholds for the difference of means were derived from the background distribution by taking the 0.5th and 99.5th percentile."
#'
#' @param eset an ExpressionSet (works with both proteins and peptides)
#' @param probs upper limit for the quantile cutoff, automatically translated to mirror the lower limit; c(1-probs, probs)
#' @param max_permutations maximum number of unique configurations used for the permuted datasets
#' @importFrom Biobase pData exprs
#' @importFrom matrixStats rowMeans2
dea_protein_background_foldchange_limits = function(eset, probs = 0.95, max_permutations = 25) {
  pd = Biobase::pData(eset)
  x = Biobase::exprs(eset)

  # find smallest group/condition
  pd_group_size = sort(table(pd$condition), decreasing = F)
  pd_ref_cond = as.numeric(names(pd_group_size)[1])
  nswap = floor(pd_group_size[1] / 2) # number of samples to change places = half of the samples from smallest group

  # all permutation of condition c (that are not the identity / unchanged)
  # this table captures the samples that should interchange labels/conditions
  pd_comb = expand.grid(apply(combn(pd$sample_id[pd$condition == pd_ref_cond], nswap), 2, paste, collapse="$$$"),
                        apply(combn(pd$sample_id[pd$condition != pd_ref_cond], nswap), 2, paste, collapse="$$$"), stringsAsFactors = F)
  # take first N permuted contrasts (or fewer if less combinations are available) to guard against computational complexity for huge data
  pd_comb = head(pd_comb, max_permutations)

  # allocate memory for all foldchanges
  fc_matrix = matrix(0.0, nrow=nrow(x), ncol=nrow(pd_comb))
  # iterate label swaps
  for(i in 1:nrow(pd_comb)) {
    # sample_id from condition 1/2 that should be swapped to condition 2/1
    sid_1 = unlist(strsplit(pd_comb[i,1], "$$$", fixed=T))
    sid_2 = unlist(strsplit(pd_comb[i,2], "$$$", fixed=T))
    # 'permuted' sample identities for each condition
    sample_id_cond1 = c(sid_2, setdiff(pd$sample_id[pd$condition==1], sid_1))
    sample_id_cond2 = c(sid_1, setdiff(pd$sample_id[pd$condition!=1], sid_2))
    # protein-level foldchange by simply taking mean value in each group
    x1 = matrixStats::rowMeans2(x[,colnames(x) %in% sample_id_cond1, drop=F])
    x2 = matrixStats::rowMeans2(x[,colnames(x) %in% sample_id_cond2, drop=F])
    fc_matrix[,i] = x1 - x2
  }

  # we cannot figure out a-symmetric foldchange thresholds in permutation data, so take the largest absolute value
  return(max(abs(quantile(fc_matrix, probs = c(1-probs, probs), na.rm = T))))

  #### some reference code for a naive random swap of sample labels. vulnerable to picking the same permutation when iterated + swapping the order within a group is meaningless
  # set.seed(1234)
  # eset_permute = eset
  # n = length(Biobase::sampleNames(eset))
  # i = sample(1:n, n)
  # pd = Biobase::pData(eset_permute)
  # pd$group = pd$group[i]
  # pd$condition = pd$condition[i]
  # Biobase::pData(eset_permute) <- pd
  # tmp = de_interface_ebayes(eset_proteins=eset, input_intensities_are_log2 = T)
  # hist(tmp$foldchange.log2)
  # tmp = de_interface_ebayes(eset_proteins=eset_permute, input_intensities_are_log2 = T)
  # quantile(tmp$foldchange.log2, probs=c(0.05, 0.995), na.rm = T)
  # hist(tmp$foldchange.log2)
}



#' convert the results from differential expression analysis from a long-format tibble to wide-format
#'
#' DEA statistics summary table in wide format: protein_id, accessions, fasta_headers, gene_symbols_or_id, <algo_de x contrast>foldchange.log2, <algo_de x contrast>pvalue, etc.
#' if there are 3 or more different DEA algorithms in the results, add a column that combines their results such that all proteins significant in 2 or more tests/algorithms are flagged
#' @param dataset your dataset. if 'de_proteins' is lacking, result is empty tibble
dea_results_to_wide = function(dataset) {
  if(!is_tibble(dataset$de_proteins) || nrow(dataset$de_proteins) == 0) {
    return(tibble())
  }

  # first, get the number of peptides used in each contrast. next, add the results from each dea algorithm in each contrast
  tib = left_join(dataset$de_proteins %>%
                    select(protein_id, contrast, peptides_used_for_dea) %>%
                    distinct(protein_id, contrast, .keep_all = T) %>%
                    pivot_wider(names_from = contrast, values_from = peptides_used_for_dea, names_prefix = "peptides_used_for_dea_") %>%
                    replace(is.na(.), 0),
                  #
                  dataset$de_proteins %>%
                    select(protein_id, algo_de, contrast, foldchange.log2, pvalue, qvalue, signif) %>%
                    pivot_wider(names_from = c(algo_de, contrast), values_from = c(foldchange.log2, pvalue, qvalue, signif)),
                  by="protein_id")

  # if there are 3 or more different DEA algorithms in the results, add a column that combines their results such that all proteins significant in 2 or more tests/algorithms are flagged
  n_algo_de = n_distinct(dataset$de_proteins$algo_de)
  if(n_algo_de > 1) {
    # from the set of significant hits, find in each contrast those protein_id that occur at least twice
    tib_signif_combined = dataset$de_proteins %>%
      filter(signif) %>%
      count(contrast, protein_id, name = "signif_count") %>%
      pivot_wider(id_cols = protein_id, names_from = contrast, values_from = signif_count, names_prefix = "signif_count_")
    # add all remaining proteins and set all NA values to zero
    tib_signif_combined = tib_signif_combined %>% right_join(tib %>% distinct(protein_id), by="protein_id") %>% replace(is.na(.), 0)

    tib = tib %>% left_join(tib_signif_combined, by = "protein_id")
  }

  return(tib)
}
