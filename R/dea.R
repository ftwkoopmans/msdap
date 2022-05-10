
#' pretty-print label for an intensity column
#' @param col todo
column_intensity_to_label = function(col) {
  ref = c("global data filter" = "intensity_all_group",
          "custom filtering and normalization" = "intensity_norm",
          "filter by group independently" = "intensity_by_group",
          "input data as-is" = "intensity")
  i = match(col, ref)
  if(!is.na(i)) {
    return(names(ref[i]))
  }

  if(grepl("^intensity_contrast:", col)) {
    return(paste0("filter by contrast;", sub("^intensity_contrast:", "", col)))
  }

  return(col)
}



#' returns named variable variable (label=column_name) indicating which type of intensity data is used
#' @param peptides todo
#' @param contr_lbl todo
get_column_intensity = function(peptides, contr_lbl = NA) {
  ref = c("filter by contrast" = paste0("intensity_", contr_lbl),
          "global data filter" = "intensity_all_group",
          "custom filtering and normalization" = "intensity_norm",
          "input data as-is" = "intensity")
  return(ref[ref %in% colnames(peptides)][1])
}



#' Differential expression analysis
#'
#' @param dataset your dataset
#' @param qval_signif threshold for significance of adjusted p-values
#' @param fc_signif threshold for significance of log2 foldchanges. Set to zero to disregard or a positive value to apply a cutoff to absolute log2 foldchanges. MS-DAP can also perform a bootstrap analyses to infer a reasonable threshold by setting this parameter to NA
#' @param algo_de algorithms for differential expression analysis. options: ebayes, deqms, msqrobsum, msempire, msqrob (to run multiple, provide an array)
#' @param algo_rollup strategy for combining peptides to proteins as used in DEA algorithms that first combine peptides to proteins and then apply statistics, like ebayes and deqms. options: maxlfq, sum. The former applies the MaxLFQ algorithm, the latter employs the 'classic' strategy of summing all peptides per protein. See further rollup_pep2prot()
#' @param output_dir_for_eset optionally, provide an output directory where the expressionset objects should be stored. Only useful if you're doing downstream data analysis that requires this data
#' @export
dea = function(dataset, qval_signif = 0.01, fc_signif = 0, algo_de = "deqms", algo_rollup = "maxlfq", output_dir_for_eset = "") {
  ### input validation
  if (length(qval_signif) != 1 || !is.finite(qval_signif) || qval_signif <= 0) {
    append_log("q-value threshold must be a single numerical value above 0 (parameter qval_signif)", type = "error")
  }

  if (length(fc_signif) != 1 || (!is.na(fc_signif) && !is.finite(fc_signif))) {
    append_log("log2 foldchange threshold must be a single numerical value or NA (parameter fc_signif)", type = "error")
  }

  if (length(algo_rollup) != 1 || ! algo_rollup %in% c("sum", "maxlfq_diann", "maxlfq")) {
    append_log("algo_rollup parameter only supports 'sum', 'maxlfq' and 'maxlfq_diann'", type = "error")
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

  # valid DEA options are those hardcoded, or pre-existing functions
  global_func = ls(envir=.GlobalEnv)
  algo_de_invalid = setdiff(algo_de, c("ebayes", "msempire", "msqrob", "msqrobsum", "deqms", global_func))
  if (length(algo_de_invalid) > 0) {
    append_log(paste("invalid options for algo_de:", paste(algo_de_invalid, collapse=", ")), type = "error")
  }

  column_contrasts = dataset_contrasts(dataset)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential expression analysis is cancelled", type = "warning")
    return(dataset)
  }


  ### data checks out
  # for computational efficiency, check whether we need protein-level data in any of the statistical models (eg; even is MaxLFQ is used as the rollup strategy, we don't want to wait for that if we are only doing msqrob)
  need_protein_rollup = any(! algo_de %in% c("msqrob","msempire"))
  # if(need_protein_rollup) {
  append_log(paste("peptide to protein rollup strategy:", algo_rollup), type = "info")
  # }

  # remove preexisting results
  dataset$de_proteins = NULL
  # init result tibble
  result_stats = tibble()

  ### iterate contrasts
  for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
    append_log(paste("differential expression analysis for", col_contr), type = "info")

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
      select(sample_id, shortname, group, condition = !!col_contr, everything()) %>%
      filter(condition != 0) %>%
      arrange(condition)

    if(!all(samples_for_contrast$sample_id %in% dataset$peptides$sample_id)) {
      append_log("all samples from this contrast must be in the peptides tibble", type = "error")
    }

    ## determine which random variables can be added to this contrast
    ranvars = ranvars_matrix = NULL
    for(v in unique(dataset$dea_random_variables)) {
      # v is a column name in dataset$samples (here we use the subset thereof for this contrast), x are the metadata values for respective column
      x = samples_for_contrast %>% pull(!!v)
      xu = unique(x) # don't re-arrange/sort !
      xi = match(x, xu)
      # because `samples_for_contrast` table was sorted by condition upstream, and condition are either 1 or 2, we can compare the random variables `x` using match(x, unique(x))
      # debug; print(cbind(samples_for_contrast$condition, match(x, unique(x)), x))
      x_aligns_with_condition = xi == samples_for_contrast$condition

      if(
        # 1) drop ranvars that align with the condition
        # is the 'random variable' not the same as the condition/contrast? criterium: at least 10% of all samples in contrast must differ (or 2)
        # May occur if user is testing many contrasts, one of which is a minor subset of the data (for which there are no differences in for instance the sample batch)
        length(xu) > 1 && sum(!x_aligns_with_condition) >= max(2, ceiling(nrow(samples_for_contrast) * .1)) &&
        # 2) drop duplicate ranvars. test the indices `xi` against ranvars from previous iterations (apply function to each column in matrix, test if all elements overlap with `xi`)
        # we must check within contrast and cannot completely check while specifying ranvars upsteam, e.g. perhaps 2 metadata properties overlap in a subset of samples that are tested in some contrast
        (length(ranvars_matrix) == 0 || !any(apply(ranvars_matrix, 2, function(col) all(col==xi))))
      ) {
        ranvars = c(ranvars, v)
        ranvars_matrix = cbind(ranvars_matrix, xi)
      }

      ## alternatively, test # differences per sample condition. proof-of-concept;
      #diff_per_condition = 0
      #for(cond in unique(samples_for_contrast$condition)) {
      #  rows = samples_for_contrast$condition == cond
      #  diff_per_condition = diff_per_condition + as.integer(n_distinct(x[rows]) > 1)
      #}
    }


    if(length(dataset$dea_random_variables) > 0 && length(ranvars) != length(dataset$dea_random_variables)) {
      append_log(paste("random variables that are _not_ applicable for current contrast due to lack of unique values (compared to sample groups/condition): ", paste(setdiff(dataset$dea_random_variables, ranvars), collapse=", ")), type = "info")
    }
    if(length(ranvars) > 0) {
      append_log(paste("random variables used in current contrast: ", paste(ranvars, collapse=", ")), type = "info")
    }


    ## convert our long-format peptide table to a peptide- and protein-level ExpressionSet
    # subset peptide tibble for current contrast
    peptides_for_contrast = dataset$peptides %>%
      select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, detect, intensity=!!as.character(col_contr_intensity), any_of("confidence")) %>%
      filter(sample_id %in% samples_for_contrast$sample_id & is.finite(intensity))
    # peptide ExpressionSet
    eset_peptides = tibble_as_eset(peptides_for_contrast, dataset$proteins, samples_for_contrast)
    # rollup peptide abundance matrix to protein-level
    m = rollup_pep2prot(peptides_for_contrast, intensity_is_log2 = T, algo_rollup = algo_rollup, return_as_matrix = T)
    # align columns with peptide-level ExpressionSet
    m = m[,match(colnames(Biobase::exprs(eset_peptides)), colnames(m)),drop=F] # use match() instead of direct key/string-based indexing because some samples may have names like 1,2,3,4 (eg; if key_sample is used for column names instead of sample_id, as we do in filter_dataset() )
    # protein ExpressionSet
    eset_proteins = protein_eset_from_data(m, eset = eset_peptides)
    rm(m)

    # eset_proteins = NULL
    # if(need_protein_rollup) {
    #   # depending on the peptide-to-protein rollup strategy, combine peptides to protein (eg; summation, maxlfq, etc.)
    #   if(algo_rollup == "maxlfq") {
    #     start_time_maxlfq = Sys.time()
    #     # maxlfq will crash on zero values. we simply threshold at 1, _assuming_ upstream data providers (eg; intensity values from user input dataset) are well above zero and any incidental value < 1 is caused by normalization of values that were already near zero
    #     peptides_for_contrast$intensity[!is.na(peptides_for_contrast$intensity) & peptides_for_contrast$intensity < 1] = 1
    #     # use the MaxLFQ implementation provided by the DIA-NN team
    #     x = diann::diann_maxlfq(peptides_for_contrast, sample.header = "sample_id", group.header="protein_id", id.header = "peptide_id", quantity.header = "intensity")
    #     x = x[, colnames(MSnbase::exprs(eset_peptides))] # align data matrix with the peptide-level ExpressionSet
    #     eset_proteins = protein_eset_from_data(x, eset = eset_peptides) # wrap the protein*sample abundance matrix in an ExpressionSet, re-using metadata from the peptide-level ExpressionSet
    #     rm(x)
    #     append_log_timestamp("peptide to protein rollup with MaxLFQ", start_time_maxlfq) # maxlfq is pretty slow, so print timer to console so the users are aware this is a time consuming step
    #   } else {
    #     eset_proteins = eset_from_peptides_to_proteins(eset_peptides, mode = algo_rollup)
    #   }
    #
    #   ## test code; optionally, normalize after non-standard rollup. very similar results when applying modebetween_protein on peptide-level data prior
    #   # if(algo_rollup == "ipl") {
    #   #   Biobase::exprs(eset_proteins) = normalize_matrix(Biobase::exprs(eset_proteins), algorithm = "vwmb", mask_sample_groups = Biobase::pData(eset_proteins)$condition)
    #   # }
    #   # if(algo_rollup == "maxlfq") {
    #   #   Biobase::exprs(eset_proteins) = normalize_matrix(Biobase::exprs(eset_proteins), algorithm = "modebetween", mask_sample_groups = Biobase::pData(eset_proteins)$condition) # post-hoc correction between groups
    #   # }
    # } else {
    #   eset_proteins = eset_from_peptides_to_proteins(eset_peptides, mode = "sum")
    # }

    # if a directory for file storage was provided, store eset in a .RData file
    if(length(output_dir_for_eset) == 1 && !is.na(output_dir_for_eset) && nchar(output_dir_for_eset)>0 && dir.exists(output_dir_for_eset)) {
      save(eset_peptides, file=path_append_and_check(output_dir_for_eset, paste0("ExpressionSet_peptides_", gsub("\\W+", " ", col_contr), ".RData")), compress = T)
      save(eset_proteins, file=path_append_and_check(output_dir_for_eset, paste0("ExpressionSet_proteins_", gsub("\\W+", " ", col_contr), ".RData")), compress = T)
    }

    contr_fc_signif = fc_signif
    if(is.na(fc_signif)) {
      contr_fc_signif = dea_protein_background_foldchange_limits(eset_proteins)
      # round the foldchange cutoff so users can get the the exact same results when they use the reported value
      contr_fc_signif = round(contr_fc_signif, digits = 3)
      append_log(sprintf("log2 foldchange threshold estimated by bootstrap analysis: %s", contr_fc_signif), type = "info")
    }

    # DE statistics for all requested algorithms
    tib = tibble()
    for(alg in algo_de) {
      err = tryCatch({
        alg_matched = F
        if (alg == "ebayes") {
          # DEBUG_esetprot <<- eset_proteins
          # DEBUG_esetpep <<- eset_peptides
          tib = bind_rows(tib, de_ebayes(eset_proteins=eset_proteins, input_intensities_are_log2 = T, random_variables = ranvars))
          alg_matched = T
        }
        if (alg == "deqms") {
          tib = bind_rows(tib, de_deqms(eset_proteins=eset_proteins, peptides=peptides_for_contrast, input_intensities_are_log2 = T, random_variables = ranvars))
          alg_matched = T
        }
        if (alg == "msempire") {
          tib = bind_rows(tib, de_msempire(eset_peptides, input_intensities_are_log2 = T)) # ! compared to the other dea algorithms, MS-EmpiRe is not a regression model so we cannot add random variables
          alg_matched = T
        }
        if (alg == "msqrob") {
          tib = bind_rows(tib, de_msqrobsum_msqrob(eset_peptides, use_peptide_model = T, input_intensities_are_log2 = T, random_variables = ranvars))
          alg_matched = T
        }
        if (alg == "msqrobsum") {
          tib = bind_rows(tib, de_msqrobsum_msqrob(eset_peptides, use_peptide_model = F, input_intensities_are_log2 = T, random_variables = ranvars))
          alg_matched = T
        }
        # for non-hardcoded functions, we call the function requested as a user parameter and pass all available data
        if(!alg_matched) {
          alg_fun = match.fun(alg)
          alg_result = alg_fun(peptides=peptides_for_contrast, samples=samples_for_contrast, eset_peptides=eset_peptides, eset_proteins=eset_proteins, input_intensities_are_log2 = T, random_variables = ranvars, dataset_name=dataset$name)
          # validation checks on expected output, to facilitate debugging/feedback for custom implementations
          if(!is_tibble(alg_result) || nrow(alg_result) == 0 || !all(c("protein_id", "pvalue", "qvalue", "foldchange.log2", "algo_de") %in% colnames(alg_result))) {
            append_log(sprintf("provided custom function for differential expression analysis '%s' must return a non-empty tibble with the columns protein_id|pvalue|qvalue|foldchange.log2|algo_de", alg), type = "error")
          }
          alg_result_name = unique(alg_result$algo_de)
          # add this criterion if method can only generate 1 result (i.e. without, 1 call can generate DEA output for my_algo_param1, my_algo_param2, etc.); length(alg_result_name) != 1 ||
          if(!is.character(alg_result_name) || nchar(alg_result_name) < 2 || alg_result_name %in% c("ebayes", "msempire", "msqrob", "msqrobsum", "combined")) {
            append_log(sprintf("provided custom function for differential expression analysis '%s' contains invalid values in algo_de column (must be a single non-empty character string uniquely indicating the name/label of your method. cannot be either of ebayes|msempire|msqrob|msqrobsum|combined)", alg), type = "error")
          }
          if(!all(!is.na(alg_result$protein_id) & alg_result$protein_id %in% peptides_for_contrast$protein_id)) {
            append_log(sprintf("provided custom function for differential expression analysis '%s' contains invalid values in protein_id column (either NA or not found in provided data structures)", alg), type = "error")
          }
          if(!all( (is.na(alg_result$pvalue) | is.numeric(alg_result$pvalue)) &
                   (is.na(alg_result$qvalue) | is.numeric(alg_result$qvalue)) &
                   (is.na(alg_result$foldchange.log2) | is.numeric(alg_result$foldchange.log2)) ))  {
            append_log(sprintf("provided custom function for differential expression analysis '%s' contains invalid values in pvalue|qvalue|foldchange.log2 columns (can only be NA or numeric)", alg), type = "error")
          }
          # finally, concatenate results
          tib = bind_rows(tib, alg_result)
        }
      }, error = function(e) e)

      if(inherits(err, "error")) {
        append_log(paste("an error occurred during the execution of DEA algorithm:", alg), type="warning")
        append_log(err$message, type="warning")
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

      result_stats = bind_rows(result_stats,
                               tib %>% left_join(prot_pep_count, by="protein_id") %>%
                                 add_column(signif = s, signif_threshold_qvalue = qval_signif, signif_threshold_log2fc = contr_fc_signif, contrast = col_contr, .after = "qvalue"))
    }
  }

  dataset$de_proteins = result_stats
  return(dataset)
}



#' estimate a threshold for 'significant' foldchanges from N permutations of sample-to-condition assignments
#'
#' note; this function hardcodes set.seed()
#'
#' Permutations of sample labels within a group are disregarded as these have no effect on the between-group foldchange, only unique combinations of swapping samples between conditions A and B are considered
#'
#' This is somewhat similar to the method described by Hafemeister and Satija at https://doi.org/10.1186/s13059-019-1874-1
#' M&M quote from this reference: "A random background distribution of mean differences was generated by randomly choosing 1000 genes and permuting the group labels. Significance thresholds for the difference of means were derived from the background distribution by taking the 0.5th and 99.5th percentile."
#'
#' @param eset an ExpressionSet (works with both proteins and peptides)
#' @param probs upper limit for the quantile cutoff, automatically translated to mirror the lower limit; c(1-probs, probs)
#' @param max_permutations maximum number of unique configurations used for the permuted datasets
#' @importFrom Biobase pData exprs
#' @importFrom matrixStats rowMeans2
dea_protein_background_foldchange_limits = function(eset, probs = 0.95, max_permutations = 100) {
  stopifnot(probs > 0.6 & probs < 1)
  stopifnot(max_permutations > 10 & max_permutations < 10000)
  pd = Biobase::pData(eset)
  x = Biobase::exprs(eset)

  ## extract sample group indentities
  # ignore samples where condition not equals 1 or 2
  samples_cond1 = pd$sample_id[pd$condition == 1]
  samples_cond2 = pd$sample_id[pd$condition == 2]
  #
  samples_cond12 = c(samples_cond1, samples_cond2)
  N1 = length(samples_cond1)
  N2 = length(samples_cond2)

  ## yields a matrix of between-group permuted indices
  index_permuted = permute_ab(c(rep(1L, N1), rep(2L, N2)), nmax = max_permutations)

  ## compute foldchange distributions from permuted sample labels
  # allocate memory for all foldchanges
  fc_matrix = matrix(0.0, nrow=nrow(x), ncol=nrow(index_permuted))
  for(i in 1:nrow(index_permuted)) { #i=1
    # extract permuted sample IDs
    samples_permuted = samples_cond12[index_permuted[i,]] # all sample IDs, re-indexed accoring to current row in permutation matrix
    sample_id_cond1 = head(samples_permuted, N1)          # first N1 samples are from group 1
    sample_id_cond2 = tail(samples_permuted, N2)          # last N2 samples are from group 2
    # print(sprintf("%d  A:%s  B:%s", i, paste(sample_id_cond1, collapse=","), paste(sample_id_cond2, collapse=",")))

    # foldchange by simply taking mean value in each group
    x1 = matrixStats::rowMeans2(x[,colnames(x) %in% sample_id_cond1, drop=F], na.rm = T)
    x2 = matrixStats::rowMeans2(x[,colnames(x) %in% sample_id_cond2, drop=F], na.rm = T)
    fc_matrix[,i] = x1 - x2
  }

  # we should not infer a-symmetric foldchange thresholds from the permutation data, so take the largest absolute value
  return(max(abs(quantile(fc_matrix, probs = c(1-probs, probs), na.rm = T))))
}



#' basic between-group permutation implementation
#'
#' minimal code & dependency-free. Not computationally efficient, but this is negligible for our use-cases
#'
#' permute_ab(as.integer(c(1,1,1, 2,2,2))) # most basic
#' permute_ab(as.integer(c(1,1,1,1, 2,2,2,2,2))) # uneven groups
#'
#' @param x integer vector of {1,2} group IDs, at least two of each
#' @param nmax maximum number of permuted sets
#' @return integer matrix where rows are permutations and columns are indices in input x
permute_ab = function(x, nmax = 100) {
  # input must be an array of only values 1 and 2, each must occur at least twice
  stopifnot(is.integer(x))
  stopifnot(x %in% 1:2)
  stopifnot(1:2 %in% x)
  stopifnot(table(x) > 1)
  stopifnot(is.numeric(nmax) & nmax > 0)

  nmax = floor(nmax) # if numeric/float, ensure this is an int

  # set random seed for reproducibility
  set.seed(1234)
  # determine group membership of each element of x
  index1 = which(x == 1)
  index2 = which(x == 2)

  ### number of elements to swap around
  # if both groups are equal size, we can get by with floor(m/2) which is slightly more efficient.
  # eg; both entry 1 and 2 have 3 observations -->> swapping 1 elements from group A to B is the same as swapping 2
  m = min(length(index1), length(index2))
  if(length(index1) == length(index2)) {
    n_swap = floor(m / 2)
  } else {
    n_swap = ceiling(m / 2)
  }

  ### naive sampling of permuted sets, overkill for very small permutations but works for any group/set sizes
  # (e.g. both groups have size=3, permutations are limited and we could just get exact solution)
  Niter = nmax * 10
  sets = matrix(0L, nrow=Niter, ncol=2*n_swap) # integer matrix, first half of the columns are sample indexes from group1
  for(i in 1:Niter) {
    # sorting is important, within-group permutations are not considered 'unique sets' in this context
    sets[i,] = c(sort(sample(index1, n_swap)), sort(sample(index2, n_swap)))
  }
  # take first N non-duplicated rows
  sets = sets[head(which(!duplicated(sets)), nmax), ]
  # debug; table(duplicated(sets)); sets[order(sets[,1],sets[,2],sets[,3],sets[,4]),]

  ### translate to permutations of x
  y = matrix(seq_along(x), nrow=nrow(sets), ncol=length(x), byrow=T)
  for(i in 1:nrow(y)) { # debug; head(y); head(sets); i=1;j=1
    for(j in 1:n_swap) {
      # swap index from group A to group B, and vice versa
      y[i,sets[i,j]] = sets[i,j+n_swap]
      y[i,sets[i,j+n_swap]] = sets[i,j]
    }
  }

  return(y)
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

  # if there are multiple DEA algorithms in the results, add a column that combines their results such that all proteins significant in 2 or more tests/algorithms are flagged
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
