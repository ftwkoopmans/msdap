
#' Plot the variance explained of sample metadata properties in the protein-intensity matrix
#'
#' precondition: the input dataset$peptides tibble must have a "intensity_all_group" column,
#' e.g. obtained by first calling the filter_dataset() function with parameter all_group=TRUE
#' (when using this function to analysis_quickstart() that is automatically taken care of)
#'
#' 1) performs peptide-to-protein rollup
#' 2) impute missing values using the missForest R package (doi: 10.1093/bioinformatics/btr597)
#' 3) call the variancePartition R package (Hoffman GE, Schadt EE, 2016, PMID:27884101)
#'
#' note; this is quite slow even on small datasets
#'
#' @param dataset your dataset. Prior to calling the this function, you must have applied "all_group" filtering and normalization using the filter_dataset() function  (see example below)
#' @param cols_metadata columns from `dataset@samples`` to be used. Set to NA or NULL to automatically infer (default)
#' @param rollup_algorithm algorithm for combining peptides to proteins as used in DEA algorithms that require a priori rollup from peptides to a protein-level abundance matrix before applying statistics (e.g. ebayes, deqms). Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#' @param quiet boolean value, passed to variancePartition package
#' @examples
#' \dontrun{
#'   dataset = filter_dataset(
#'     dataset, filter_min_detect = 0, filter_min_quant = 3,
#'     by_group = FALSE, by_contrast = FALSE, all_group = TRUE,
#'     norm_algorithm = c("vwmb", "modebetween_protein")
#'   )
#'   tmp = plot_variance_explained(dataset, cols_metadata = NA)
#'   print(tmp$p_ve_violin)
#'   print(tmp$tbl_ve)
#' }
#'
#' @importFrom variancePartition fitExtractVarPartModel plotVarPart
#' @importFrom missForest missForest
#' @importFrom BiocParallel SerialParam
#' @export
plot_variance_explained = function(dataset, cols_metadata = NULL, rollup_algorithm = "maxlfq", quiet = TRUE) {
  start_time = Sys.time()
  if(! "intensity_all_group" %in% colnames(dataset$peptides)) {
    append_log("'all group' filtering has to be applied before running this function (e.g. filter_dataset(..., all_group=TRUE) )", type = "error")
  }

  ## sample metadata columns
  cols_metadata = setdiff(na.omit(cols_metadata), "") # from input, disregard NA and empty string
  cols_metadata_auto = user_provided_metadata(dataset$samples)

  if(length(cols_metadata) == 0) { # if none provided, infer from sample metadata table
    cols_metadata = cols_metadata_auto
  } else { # if provided, check that valid options are used (i.e. catch non-existing columns, or auto-generated metadata like sample_index)
    if(! all(cols_metadata %in% cols_metadata_auto)) {
      append_log(paste0("invalid metadata columns provided @ plot_variance_explained(). Valid options: ", paste(cols_metadata_auto, collapse = ","), ". Invalid input provided: ", paste(setdiff(cols_metadata, cols_metadata_auto)), collapse = ","), type = "error")
    }
  }

  if(length(cols_metadata) == 0) {
    append_log("no sample metadata available for variance-explained analysis", type = "warning")
    return() # no metadata, no plot
  }


  ## protein-level matrix
  protein_matrix = rollup_pep2prot(
    dataset$peptides %>%
      filter(!is.na(intensity_all_group)) %>%
      select(sample_id, protein_id, peptide_id, intensity=intensity_all_group),
    intensity_is_log2 = T,
    rollup_algorithm = rollup_algorithm,
    return_as_matrix = T)


  ## extract sample metadata
  s = dataset$samples %>%
    # only samples that remain after the filter @ intensity_all_group
    filter(sample_id %in% colnames(protein_matrix)) %>%
    select(sample_id, !!cols_metadata)
  # for the variancePartition package we need to rename the sample ID column as 'ID'. So if user metadata has that name as well, rename it
  if("ID" %in% colnames(s)) {
    s = s %>% rename(ID_ = ID)
  }
  s = s %>% rename(ID = sample_id)

  ## convert character columns to factors
  s_types = sapply(s, class) # don't use typeof(), it returns 'integer' for factors
  cols_fix = setdiff(names(s_types)[s_types %in% c("character", "logical")], "ID")
  cols_na_repaired = NULL
  for(n in cols_fix) {
    tmp = s %>% select(any_of(n)) %>% pull() %>% as.character() # grab values from column and convert to string
    if(n_distinct(na.omit(tmp)) < 2) {
      append_log(paste("sample metadata column that only contains a single value:", n), type = "error")
    }
    if(anyNA(tmp)) {
      cols_na_repaired = c(cols_na_repaired, n)
      tmp[is.na(tmp)] = "NA" # replace NA with a string value
    }
    s[,n] = factor(tmp, levels = unique(tmp)) # factor levels in same order as input table, i.e. don't auto-sort factor levels
  }

  ## numeric columns cannot be NA !
  cols_numeric = setdiff(names(s_types)[s_types %in% c("numeric", "integer")], "ID")
  for(n in cols_numeric) {
    tmp = s %>% select(any_of(n)) %>% pull()
    if(n_distinct(na.omit(tmp)) < 2) {
      append_log(paste("sample metadata column that only contains a single value:", n), type = "error")
    }
    if(any(!is.finite(tmp))) {
      cols_na_repaired = c(cols_na_repaired, n)
      tmp[!is.finite(tmp)] = mean(tmp, na.rm = T)
      s[,n] = tmp
    }
  }

  if(length(cols_na_repaired) > 0) {
    append_log(paste0("plot_variance_explained: replaced NA values in these sample metadata columns; ", paste(cols_na_repaired, collapse = ","), "\nfor numeric columns the mean value was imputed, for other column types the missing values were replaced by the string 'NA' (e.g. a boolean sample metadata column with missing values would become a factor with levels 'TRUE'/'FALSE'/'NA')"), type = "warning")
  }



  append_log("computing variance-explained...", type = "info")

  ## DDA data may have missing values and unfortunately the variancePartition package cannot deal with this
  # assuming peptide filters were applied upstream, the % missingness is low and imputation has relatively low impact
  # we here use Random Forest imputation @ missForest R package, which should have minimal impact on this analysis (again, provided there are relatively few missing values)
  # Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
  # reference for chosing this imputation approach: Jin L, Bi Y, Hu C, Qu J, Shen S, Wang X, Tian Y. A comparative study of evaluating missing value imputation methods in label-free proteomics. Sci Rep. 2021 Jan 19;11(1):1760. doi: 10.1038/s41598-021-81279-4. PMID: 33469060
  if(anyNA(protein_matrix)) {
    capture.output( protein_matrix <- missForest::missForest(protein_matrix, verbose = F)$ximp ) # note, must use <- notation within capture.output
  }


  ## apply variancePartition package
  # construct the formula. note that the "ID" column should be excluded
  formula_string__factors = sprintf("(1 | %s)", setdiff(colnames(s), c("ID", cols_numeric)))
  formula_string = paste("~", paste(c(formula_string__factors, cols_numeric), collapse=" + "))
  append_log(paste("formula for variancePartition;", formula_string), type = "info")
  # formula_string = paste("~", paste(sprintf("(1 | %s)", setdiff(colnames(s), "ID")), collapse=" + "))

  vp = tryCatch({
    ## alas, parallelism is very buggy for R 3.6 on Windows. workarounds below
    # tell foreach to use a 1 thread sequential worker, so the internal function `variancePartition:::.isDisconnected()` , which runs foreach::foreach(i = seq_len(2)) %dopar% { i } as QC check, doesn't yield an error if any upstream-pre-registered SNOW cluster bugs out
    foreach::registerDoSEQ()
    suppressWarnings(suppressMessages(variancePartition::fitExtractVarPartModel(
      protein_matrix, formula = eval(parse(text=formula_string)), data = s, quiet = quiet, showWarnings = F,
      # hardcoded low-N parallelization to counter known bug; http://bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/FAQ.html#errors-with-biocparallel-multithreading-backend
      # instead of SNOW as suggested in the manual / link above, we disabled parallelism altogether as it is unreliable across some R 3.6.3 Windows systems we tested  (roll this out to wider public = ask for trouble)
      BPPARAM = BiocParallel::SerialParam()
    )))
  },
  error = function(err) {
    append_log(paste("failed to run variance-explained analysis. Often this is due to rank deficiency while fitting the linear models, there is insufficient information contained in your data to estimate the model created from all metadata parameters;", formula_string, " -->> try again by specifying a few non-redundant properties that represent categorical variables (e.g. set the 'cols_metadata' parameter to only:  c('group',  'batch') )"), type = "warning")
    return(NULL)
  })

  if(length(vp) > 0) {
    append_log_timestamp("variance-explained", start_time)

    ## collect plot and summary stats (as plot/table)
    vp_sort = variancePartition::sortCols(variancePartition::sortCols(vp, FUN = mean), fun = stats::median)
    p_ve_violin = variancePartition::plotVarPart(vp_sort)
    tbl_ve = as_tibble(as.matrix(vp_sort)) %>%
      summarize_all(.funs = function(x) {bp=grDevices::boxplot.stats(x, do.conf=FALSE, do.out=FALSE)$stats; m=mean(x,na.rm=T); c(bp[5:4], m, bp[3:1])} ) %>% mutate_all(function(x) sprintf("%.1f", x * 100)) %>%
      add_column(statistic = c("boxplot upper whisker", "boxplot upper hinge", "mean", "median", "boxplot lower hinge", "boxplot lower whisker"), .before = 1)

    return(list(p_ve_violin=p_ve_violin, tbl_ve = tbl_ve, ve_data = vp_sort))
  }
}

