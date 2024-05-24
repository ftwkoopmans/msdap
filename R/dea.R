
#' returns all DEA functions integrated with MS-DAP
#'
#' @description
#'
#' ## available DEA functions
#' **ebayes**: wrapper for the eBayes function from the limma package (PMID:25605792) <https://bioconductor.org/packages/release/bioc/html/limma.html>. The eBayes function applies moderated t-statistics to each row of a given data matrix.
#' It was originally developed for the analysis of RNA-sequencing data but can also be applied to a proteomics `protein*sample` data matrix. Doesn't work on peptide-level data because limma eBayes returns t-statistics per row in the input matrix,
#' so using peptide-level data would yield statistics per peptide and translating those to protein-level statistics is not straight forward.
#' Thus, with MS-DAP we first perform peptide-to-protein rollup (e.g. using MaxLFQ algorithm) and then apply the limma eBayes function to the protein-level data matrix.
#'
#' This is in line with typical usage of moderated t-statistics in proteomics (e.g. analogous to Perseus, where a moderated t-test is applied to protein data matrix). This method will take provided covariates into account (if any). Implemented in function; \code{de_ebayes}
#'
#' **deqms**: wrapper for the DEqMS package (PMID:32205417) <https://github.com/yafeng/DEqMS>, which is a proteomics-focussed extension of the limma eBayes approach that also weighs the number of peptides observed per protein to adjust protein statistics. MS-DAP will apply this function to the protein-level data matrix. This method will take provided covariates into account (if any). Implemented in function; \code{de_deqms}
#'
#' **msempire**: wrapper for the msEmpiRe package (PMID:31235637) <https://github.com/zimmerlab/MS-EmpiRe>. This is a peptide-level DEA algorithm. Note that this method cannot deal with covariates! Implemented in function; \code{de_msempire}
#'
#' **msqrob**: implementation of the MSqRob package, with minor tweak for (situationally) faster computation (PMID:26566788) <https://github.com/statOmics/msqrob>. This is a peptide-level DEA algorithm. This method will take provided covariates into account (if any). Implemented in function; \code{de_msqrobsum_msqrob}
#'
#' **msqrob_fc**: msqrob estimated protein pvalues, but with log2fc values computed by eBayes. This avoids the excessive shrinkage of protein log2fc values by MSqRob, but does use the MSqRob robust peptide regression model for protein p-value estimates (which outperforms eBayes/DEqMS in ROC analyses on most benchmark datasets, so this method is intended as a "best of both worlds")
#'
#' **msqrobsum**: implementation of the MSqRob package (which also features MSqRobSum), with minor tweak for (situationally) faster computation (PMID:32321741) <https://github.com/statOmics/msqrob>. This is a hybrid peptide&protein-level DEA algorithm that takes peptide-level data as input; it first performs peptide-to-protein rollup, then applies statistics to this protein-level data matrix. This method will take provided covariates into account (if any). Implemented in function; \code{de_msqrobsum_msqrob}
#'
#' @export
dea_algorithms = function() {
  return(c("ebayes", "deqms", "msempire", "msqrobsum", "msqrob")) # , "msqrob_fc"
}



#' reserved DEA function names
#'
#' The returned function names are reserved; custom DEA functions cannot have any of these names
#' @export
dea_algorithms_reserved = function() {
  return(c(dea_algorithms(), "combined"))
}



#' pretty-print label for an intensity column
#'
#' @param col column name in dataset$peptides
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



#' prioritized selection of intensity column in a peptides tibble
#'
#' 1) if a column with intensity values for a specified contrasts exists, return it
#' 2) all_group filtering in "intensity_all_group" column
#' 3) custom filtering in "intensity_norm" column
#' 4) raw data in "intensity" column
#'
#' returns a named array of length 1 where the name is a "pretty print label" and value is the column name
#'
#' @param peptides peptides tibble, e.g. dataset$peptides
#' @param contr_lbl optionally a string describing a contrast, e.g. "contrast: WT vs KO"
get_column_intensity = function(peptides, contr_lbl = NA) {
  ref = c("filter by contrast" = paste0("intensity_", contr_lbl),
          "global data filter" = "intensity_all_group",
          "custom filtering and normalization" = "intensity_norm",
          "input data as-is" = "intensity")
  return(ref[ref %in% colnames(peptides)][1])
}



#' return sample groups for a given contrast (column name in dataset$samples)
#'
#' @param dataset your dataset
#' @param contr_lbl a string describing a contrast = column in dataset$samples, e.g. "contrast: WT vs KO"
#' @returns list of unique sample groups, first element is left-hand side of the contrast and second element the right-hand side
contrast_to_samplegroups = function(dataset, contr_lbl) {
  stopifnot(contr_lbl %in% colnames(dataset$samples))
  x = dataset$samples %>% select(group, contrast = !!contr_lbl)
  return(list(unique(x$group[x$contrast == 1]), unique(x$group[x$contrast == 2])))
}



#' Differential expression analysis
#'
#' @param dataset your dataset
#' @param qval_signif threshold for significance of q-values
#' @param fc_signif threshold for significance of log2 foldchanges. Set to zero or NA to disregard, or a positive value to apply a cutoff to absolute log2 foldchanges. MS-DAP can also perform a bootstrap analyses to infer a reasonable threshold by setting this parameter to NA
#' @param dea_algorithm algorithm for differential expression analysis (provide an array of strings to run multiple, in parallel). Refer to \code{\link{dea_algorithms}} function documentation for available options and a brief description of each. To use a custom DEA function, provide the respective R function name as a string (see GitHub documentation on custom DEA functions for more details)
#' @param rollup_algorithm algorithm for combining peptides to proteins as used in DEA algorithms that require a priori rollup from peptides to a protein-level abundance matrix before applying statistics (e.g. ebayes, deqms). Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#' @param output_dir_for_eset optionally, provide an output directory where the expressionset objects should be stored. Only useful if you're doing downstream data analysis that requires this data
#' @return long-format tibble with results for each DEA algorithm requested via `dea_algorithm` parameter. Note that a MS-DAP contrast for "A vs B" returns foldchanges for B/A. For example, for the contrast "control vs disease" a positive log2 foldchange implies protein abundances are higher in the "disease" sample group. The column `signif` contains a boolean flag indicating significant hits according to both the user defined q-value threshold (parameter `qval_signif`, also in result table column `signif_threshold_qvalue`) and optional foldchange threshold (parameter `fc_signif`, also in result table column `signif_threshold_log2fc`)
#' @seealso `dea_algorithms()` for available DEA algorithms and documentation.
#' @export
dea = function(dataset, qval_signif = 0.01, fc_signif = 0, dea_algorithm = "deqms", rollup_algorithm = "maxlfq", output_dir_for_eset = "") {
  ### input validation
  if (length(qval_signif) != 1 || !is.finite(qval_signif) || qval_signif <= 0) {
    append_log("q-value threshold must be a single numerical value above 0 (parameter qval_signif)", type = "error")
  }

  if (length(fc_signif) != 1 || (!is.na(fc_signif) && !is.finite(fc_signif))) {
    append_log("log2 foldchange threshold must be a single numerical value or NA (parameter fc_signif)", type = "error")
  }

  if (length(rollup_algorithm) != 1 || ! rollup_algorithm %in% c("sum", "maxlfq_diann", "maxlfq", "tukey_median")) {
    append_log("rollup_algorithm parameter only supports 'sum', 'maxlfq', 'maxlfq_diann' and 'tukey_median'", type = "error")
  }


  # trigger automatic estimation of FC by either NA or a negative value (to convenience the user, no functional implications)
  if(!is.na(fc_signif) && fc_signif < 0) {
    fc_signif = NA
  }

  dea_algorithm = unique(dea_algorithm)
  if (length(dea_algorithm) == 0) {
    append_log("no algorithms have been defined (parameter dea_algorithm), differential expression analysis is cancelled", type = "warning")
    return(dataset)
  }

  # valid DEA options are those hardcoded, or pre-existing functions
  global_func = ls(envir=.GlobalEnv)
  dea_algorithm_invalid = setdiff(dea_algorithm, c(dea_algorithms(), global_func))
  if (length(dea_algorithm_invalid) > 0) {
    append_log(paste("invalid options for dea_algorithm:", paste(dea_algorithm_invalid, collapse=", ")), type = "error")
  }

  column_contrasts = dataset_contrasts(dataset)
  if (length(column_contrasts) == 0) {
    append_log("no contrasts have been defined, differential expression analysis is cancelled", type = "warning")
    return(dataset)
  }

  # remove preexisting results
  dataset$de_proteins = NULL
  # init result tibble
  result_stats = tibble()

  ### iterate contrasts
  for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
    append_log(paste("differential expression analysis for", col_contr), type = "info")

    # returns named variable (label=column_name) indicating which type of intensity data is used
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
      select(sample_id, shortname, group, condition = !!col_contr, tidyselect::everything()) %>%
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
      select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, detect, intensity=!!as.character(col_contr_intensity), any_of(c("confidence","charge"))) %>%
      filter(sample_id %in% samples_for_contrast$sample_id & is.finite(intensity))
    # peptide ExpressionSet
    eset_peptides = tibble_as_eset(peptides_for_contrast, dataset$proteins, samples_for_contrast)
    # rollup peptide abundance matrix to protein-level
    m = rollup_pep2prot(peptides_for_contrast, intensity_is_log2 = T, rollup_algorithm = rollup_algorithm, return_as_matrix = T)
    # align columns with peptide-level ExpressionSet
    m = m[,match(colnames(Biobase::exprs(eset_peptides)), colnames(m)),drop=F] # use match() instead of direct key/string-based indexing because some samples may have names like 1,2,3,4 (eg; if key_sample is used for column names instead of sample_id, as we do in filter_dataset() )
    # protein ExpressionSet
    eset_proteins = protein_eset_from_data(m, eset = eset_peptides)
    # count the number of unique peptides per protein and add it to the ExpressionSet protein metadata data.frame ('fData')
    prot_pep_count = peptides_for_contrast %>% distinct(protein_id, peptide_id) %>% count(protein_id)
    Biobase::fData(eset_proteins) = Biobase::fData(eset_proteins) %>% mutate(npep = prot_pep_count$n[match(protein_id, prot_pep_count$protein_id)])
    # cleanup
    rm(m, prot_pep_count)

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
    for(alg in dea_algorithm) {
      err = tryCatch({
        ### call DEA function
        alg_result = NULL
        alg_plugin = FALSE
        if(alg == "ebayes") {
          alg_result = de_ebayes(eset_proteins=eset_proteins, input_intensities_are_log2 = T, random_variables = ranvars)
        } else if(alg == "deqms") {
          alg_result = de_deqms(eset_proteins=eset_proteins, input_intensities_are_log2 = T, random_variables = ranvars)
        } else if(alg == "msempire") {
          # ! compared to the other dea algorithms, MS-EmpiRe is not a regression model so we cannot add random variables
          alg_result = de_msempire(eset_peptides, input_intensities_are_log2 = T)
        } else if(alg == "msqrob") {
          alg_result = de_msqrobsum_msqrob(eset_peptides, eset_proteins = NULL, use_peptide_model = T, input_intensities_are_log2 = T, random_variables = ranvars, log2fc_without_shrinkage = FALSE)
        # } else if(alg == "msqrob_fc") {
        #   alg_result = de_msqrobsum_msqrob(eset_peptides, eset_proteins = eset_proteins, use_peptide_model = T, input_intensities_are_log2 = T, random_variables = ranvars, log2fc_without_shrinkage = TRUE)
        } else if(alg == "msqrobsum") {
          alg_result = de_msqrobsum_msqrob(eset_peptides, eset_proteins = NULL, use_peptide_model = F, input_intensities_are_log2 = T, random_variables = ranvars, log2fc_without_shrinkage = FALSE)
        } else {
          # for non-hardcoded functions, we call the function requested as a user parameter and pass all available data
          alg_fun = match.fun(alg) # throws error if not found
          alg_result = alg_fun(peptides=peptides_for_contrast, samples=samples_for_contrast, eset_peptides=eset_peptides, eset_proteins=eset_proteins, input_intensities_are_log2 = T, random_variables = ranvars, dataset_name=dataset$name)
          alg_plugin = TRUE
        }

        ### validate results from DEA
        # non-empty result table with required columns
        if(!is_tibble(alg_result) || nrow(alg_result) == 0 || !all(c("protein_id", "pvalue", "qvalue", "foldchange.log2", "dea_algorithm") %in% colnames(alg_result))) {
          append_log(sprintf("DEA function '%s' must return a non-empty tibble with the columns protein_id|pvalue|qvalue|foldchange.log2|dea_algorithm", alg), type = "error")
        }
        # dea_algorithm ID must be a valid string and not be a reserved method name. Note that a method is allowed to return multiple (e.g. 1 wrapper function that returns stats for both algorithm 1 and algorithm 2)
        alg_result_name = unique(alg_result$dea_algorithm)
        if(length(alg_result_name) == 0 || any(!is.character(alg_result_name) | !grepl("^[0-9a-zA-Z_-]{2,}$", alg_result_name))) {
          append_log(sprintf("DEA function '%s' results contain invalid values in dea_algorithm column; must be a character string with only letters/numbers/underscores (length >= 2), uniquely indicating the name of your method (to be used in plots and output tables)", alg), type = "error")
        }
        if(alg_plugin && any(alg_result_name %in% dea_algorithms_reserved())) {
          append_log(sprintf("DEA function '%s' results contain invalid values in dea_algorithm column; cannot be either of these reserved keywords: %s", alg, paste(dea_algorithms_reserved(), collapse = ",")), type = "error")
        }
        # if, within this contrast, we've already stored results for this DEA algorithm... user passed multiple custom functions that return same dea_algorithm ID
        if(nrow(tib) > 0 && any(alg_result_name %in% tib$dea_algorithm)) {
          append_log(sprintf("DEA function '%s' results contain invalid values in dea_algorithm column; results for this dea_algorithm ID have already been stored previously. Did you provide multiple DEA functions that return the same dea_algorithm value?", alg), type = "error")
        }
        if(any(is.na(alg_result$protein_id) | !is.character(alg_result$protein_id) | ! alg_result$protein_id %in% unique(peptides_for_contrast$protein_id))) {
          append_log(sprintf("DEA function '%s' results contain invalid values in protein_id column (NA not allowed, must be characters not factors, returned protein_id values must all be present in peptide tibble provided as parameter to DEA function)", alg), type = "error")
        }
        if(!all( (is.na(alg_result$pvalue) | is.numeric(alg_result$pvalue)) & !is.infinite(alg_result$pvalue) &
                 (is.na(alg_result$qvalue) | is.numeric(alg_result$qvalue)) & !is.infinite(alg_result$qvalue) &
                 (is.na(alg_result$foldchange.log2) | is.numeric(alg_result$foldchange.log2)) &!is.infinite(alg_result$foldchange.log2) ))  {
          append_log(sprintf("DEA function '%s' results contain invalid values in pvalue|qvalue|foldchange.log2 columns (can only be NA or numeric, no Infinite allowed)", alg), type = "error")
        }
        # 1 dea_algorithm can yield only 1 result per protein_id. So test if # unique combinations of both is the same as input table length
        if(alg_result %>% distinct(dea_algorithm, protein_id) %>% nrow() != nrow(alg_result)) { # faster than `anyDuplicated(alg_result %>% select(dea_algorithm, protein_id)) != 0`
          append_log(sprintf("DEA function '%s' results contain duplicate values in protein_id column (within same dea_algorithm)", alg), type = "error")
        }

        ### finally, concatenate results
        tib = bind_rows(tib, alg_result)

      }, error = function(e) e)

      # error handling; a DEA function returned an error, show to user using our logging system
      if(inherits(err, "error")) {
        append_log(sprintf("an error occurred in %s during the execution of DEA function '%s'", col_contr, alg), type="warning")
        append_log_error(err, type="warning")
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

      result_stats = bind_rows(
        result_stats,
        tib %>%
          left_join(prot_pep_count, by="protein_id") %>%
          add_column(signif = s, signif_threshold_qvalue = qval_signif, signif_threshold_log2fc = contr_fc_signif, contrast = col_contr, .after = "qvalue")
      )
    }
  }

  if(is.data.frame(result_stats)) {
    result_stats = result_stats %>% mutate_all(unname)
  }

  dataset$de_proteins = result_stats
  return(dataset)
}



#' estimate a threshold for 'significant' foldchanges from N permutations of sample-to-condition assignments
#'
#' @description
#' note; this function hardcodes set.seed()
#'
#' Permutations of sample labels within a group are disregarded as these have no effect on the between-group foldchange, only unique combinations of swapping samples between conditions A and B are considered
#'
#' This is somewhat similar to the method described by Hafemeister and Satija at https://doi.org/10.1186/s13059-019-1874-1
#' M&M quote from this reference: "A random background distribution of mean differences was generated by randomly choosing 1000 genes and permuting the group labels. Significance thresholds for the difference of means were derived from the background distribution by taking the 0.5th and 99.5th percentile."
#' @examples \dontrun{
#' m = matrix(runif(600), nrow=100, ncol=6, dimnames=list(NULL, LETTERS[1:6]))
#' dea_protein_background_foldchange_limits__v2(
#'   m,
#'   samples_cond1 = c("A","B","C"), samples_cond2 = c("D","E","F"),
#'   probs = 0.95, max_permutations = 100
#' )
#' }
#' @param x a numeric matrix or a Biobase::ExpressionSet(). If the latter is provided, the ExpressionSet must contain a column named "condition" in phenotypic data and it must contain numeric values (which will be casted to integers), where 1 indicates condition 1 and 2 indicates condition 2
#' @param samples_cond1 if x is a matrix, this must be a character vector indicating samples from condition 1; these must be column names in x
#' @param samples_cond2 analogous to samples_cond1, but for condition 2
#' @param probs upper limit for the quantile cutoff, automatically translated to mirror the lower limit; c(1-probs, probs)
#' @param max_permutations maximum number of unique configurations used for the permuted datasets
#' @importFrom Biobase pData exprs
#' @importFrom matrixStats rowMeans2
#' @export
dea_protein_background_foldchange_limits = function(x, samples_cond1 = NULL, samples_cond2 = NULL, probs = 0.95, max_permutations = 100) {
  x_ismatrix = is.matrix(x) && is.numeric(x) && length(samples_cond1) > 0 && length(samples_cond2) > 0 && all(samples_cond1 %in% colnames(x)) && all(samples_cond2 %in% colnames(x))
  x_iseset = "ExpressionSet" %in% class(x) && "condition" %in% colnames(Biobase::pData(x)) && all(c(1L,2L) %in% as.integer(Biobase::pData(x)$condition))
  if(!x_ismatrix && !x_iseset) {
    stop("dea_protein_background_foldchange_limits() requires a valid ExpressionSet or a numeric matrix together with column-name-to-condition mappings")
  }
  stopifnot(probs > 0.6 & probs < 1)
  stopifnot(max_permutations > 10 & max_permutations < 10000)

  if(x_iseset) {
    pd = Biobase::pData(x) # extract phenotype data
    x = Biobase::exprs(x) # overwrite x with expression matrix from the eset
    # extract sample group indentities, ignoring samples where condition not equals 1 or 2
    samples_cond1 = pd$sample_id[as.integer(pd$condition) == 1L]
    samples_cond2 = pd$sample_id[as.integer(pd$condition) == 2L]

    rm(pd)
  }

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
#' DEA statistics table `dataset$de_proteins` presented in wide format;
#' protein_id, accessions, fasta_headers, gene_symbols_or_id
#' <dea_algorithm x contrast>foldchange.log2, <dea_algorithm x contrast>pvalue, <dea_algorithm x contrast>qvalue (pvalue adjusted for multiple testing), etc.
#'
#' if there are 2 or more different DEA algorithms in the results, add a column that combines their results such that all proteins significant in 2 or more tests/algorithms are flagged
#'
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
                    select(protein_id, dea_algorithm, contrast, foldchange.log2, pvalue, qvalue, signif) %>%
                    pivot_wider(names_from = c(dea_algorithm, contrast), values_from = c(foldchange.log2, pvalue, qvalue, signif)),
                  by="protein_id")

  # if there are multiple DEA algorithms in the results, add a column that combines their results such that all proteins significant in 2 or more tests/algorithms are flagged
  n_dea_algorithm = n_distinct(dataset$de_proteins$dea_algorithm)
  if(n_dea_algorithm > 1) {
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
