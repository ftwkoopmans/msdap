
#' Plot the variance explained of sample metadata properties in the protein-intensity matrix
#'
#' basically a wrapper around the variancePartition R package (Hoffman GE, Schadt EE, 2016, PMID:27884101)
#'
#' note; this is quite slow even on small datasets
#'
#' @param dataset your dataset
#' @param cols_metadata columns from dataset@samples to be used. leave empty to automatically infer (default)
#' @importFrom variancePartition fitExtractVarPartModel plotVarPart
#' @importFrom missForest missForest
#' @importFrom BiocParallel SerialParam
#' @importFrom ggpubr ggarrange
plot_variance_explained = function(dataset, cols_metadata = NULL) {
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

  append_log("computing variance-explained...", type = "info")

  ## protein-level matrix
  protein_matrix = msdap::rollup_pep2prot(
    dataset$peptides %>%
      filter(!is.na(intensity_all_group)) %>%
      select(sample_id, protein_id, peptide_id, intensity=intensity_all_group),
    intensity_is_log2 = T,
    algo_rollup = "maxlfq",
    return_as_matrix = T)

  ## DDA data may have missing values and unfortunately the variancePartition package cannot deal with this
  # assuming peptide filters were applied upstream, the % missingness is low and imputation has relatively low impact
  # we here use Random Forest imputation @ missForest R package, which should have minimal impact on this analysis (again, provided there are relatively few missing values)
  # Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
  # reference for chosing this imputation approach: Jin L, Bi Y, Hu C, Qu J, Shen S, Wang X, Tian Y. A comparative study of evaluating missing value imputation methods in label-free proteomics. Sci Rep. 2021 Jan 19;11(1):1760. doi: 10.1038/s41598-021-81279-4. PMID: 33469060
  if(any(is.na(protein_matrix))) {
    capture.output( protein_matrix <- missForest::missForest(protein_matrix, verbose = F)$ximp ) # note, must use <- notation within capture.output
  }


  ## extract sample metadata
  s = dataset$samples %>%
    # only samples that remain after the filter @ intensity_all_group
    filter(sample_id %in% colnames(protein_matrix)) %>%
    select(sample_id, !!cols_metadata)
  # for the variancePartition package we need to rename the sample ID column as 'ID'. So if user metadata has that name as well, rename it
  if("ID" %in% names(s)) {
    s = s %>% rename(ID_ = ID)
  }
  s = s %>% rename(ID = sample_id)

  # convert to character factors
  metadata = s %>%
    mutate_all(as.character) %>%
    replace(is.na(.), "NA") %>%
    mutate_all(as.factor)


  ## apply variancePartition package
  formula_string = paste("~", paste(sprintf("(1 | %s)", setdiff(colnames(metadata), "ID")), collapse=" + "))
  vp = tryCatch({
    ## alas, parallelism is very buggy for R 3.6 on Windows. workarounds below
    # tell foreach to use a 1 thread sequential worker, so the internal function `variancePartition:::.isDisconnected()` , which runs foreach::foreach(i = seq_len(2)) %dopar% { i } as QC check, doesn't yield an error if any upstream-pre-registered SNOW cluster bugs out
    foreach::registerDoSEQ()
    suppressWarnings(suppressMessages(variancePartition::fitExtractVarPartModel(protein_matrix, formula = eval(parse(text=formula_string)), data = metadata, quiet = T, showWarnings = F,
                                                                                      # hardcoded low-N parallelization to counter known bug; http://bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/FAQ.html#errors-with-biocparallel-multithreading-backend
                                                                                      # instead of SNOW as suggested in the manual / link above, we disabled parallelism altogether as it is unreliable across some R 3.6.3 Windows systems we tested  (roll this out to wider public = ask for trouble)
                                                                                      BPPARAM = BiocParallel::SerialParam()) ))
  },
  error = function(err) {
    append_log(paste("failed to run variance-explained analysis. Often this is due to rank deficiency while fitting the linear models, there is insufficient information contained in your data to estimate the model created from all metadata parameters;", formula_string, " -->> try again by hardcoding a few non-redundant properties that represent categorical variables (e.g. set the parameter to:  c('group', 'gel', 'batch') )"), type = "warning")
    return(NULL)
  })

  if(length(vp) > 0) {
    append_log_timestamp("variance-explained", start_time)

    ## collect plot and summary stats (as plot/table)
    vp_sort = variancePartition::sortCols(variancePartition::sortCols(vp, FUN = mean), fun = median)
    p_ve_violin = variancePartition::plotVarPart(vp_sort, )
    tbl_ve = as_tibble(as.matrix(vp_sort)) %>%
      summarize_all(.funs = function(x) {bp=boxplot.stats(x, do.conf=FALSE, do.out=FALSE)$stats; m=mean(x,na.rm=T); c(bp[5:4], m, bp[3:1])} ) %>% mutate_all(function(x) sprintf("%.1f", x * 100)) %>%
      add_column(statistic = c("boxplot upper whisker", "boxplot upper hinge", "mean", "median", "boxplot lower hinge", "boxplot lower whisker"), .before = 1)

    # p_ve_table = ggpubr::ggtexttable(tbl_ve, rows = NULL)
    # p = ggpubr::ggarrange(p_ve_violin, p_ve_table, ncol = 1, nrow = 2, heights = c(3,1) )

    return(list(p_ve_violin=p_ve_violin, tbl_ve = tbl_ve))
  }
}

