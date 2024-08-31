

#' Input validation for eBayes/DEqMS/MS-EmpiRe functions
#'
#' @param eset protein/peptide log2 intensity matrix stored as a Biobase ExpressionSet. Must describe protein_id and sample_id in metadata
#' @param model_matrix a `stats::model.matrix()` result that is supplied to `limma::lmFit()`
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
validate_eset_fitdata = function(eset, model_matrix, model_matrix_result_prop) {
  if(length(eset) != 1 || !"ExpressionSet" %in% class(eset)) {
    append_log("parameter eset must be an ExpressionSet", type = "error")
  }

  f = Biobase::fData(eset)
  if(!is.data.frame(f) || !"protein_id" %in% colnames(f) || anyNA(f$protein_id)) {
    append_log('parameter eset must be an ExpressionSet with non-NA "protein_id" as a featuredata (i.e. Biobase::fData(eset)$protein_id)', type = "error")
  }

  p = Biobase::pData(eset)
  if(!is.data.frame(p) || !"sample_id" %in% colnames(p) || anyNA(p$sample_id)) {
    append_log('parameter eset must be an ExpressionSet with non-NA "sample_id" as a phenotypic data  (i.e. Biobase::pData(eset)$sample_id)', type = "error")
  }

  if(!is.matrix(model_matrix)) {
    append_log("parameter model_matrix must be a matrix", type = "error")
  }

  if(nrow(model_matrix) != nrow(p)) {
    append_log("parameters eset and model_matrix must describe the same number of samples (i.e. number of columns/samples in ExpressionSet must match the number of rows in model_matrix)", type = "error")
  }

  if(length(model_matrix_result_prop) != 1 || !is.character(model_matrix_result_prop) ||
     model_matrix_result_prop == "" || !model_matrix_result_prop %in% colnames(model_matrix)) {
    append_log("parameter model_matrix_result_prop must be a single string, matching any of the column names in the model_matrix parameter", type = "error")
  }
}



#' helper function to double-check we don't submit an ExpressionSet with too few values to DEA
#'
#' @param eset a Biobase ExpressionSet
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
#' @param is_peptide boolean indicating peptide-level (FALSE = protein-level)
#' @param method_name label for printing error messages
subset_protein_eset_for_dea = function(eset, model_matrix_result_prop, is_peptide, method_name) {
  x = Biobase::exprs(eset)
  ptmp = Biobase::pData(eset)
  stopifnot("subset_protein_eset_for_dea invalid eset * model_matrix_result_prop combination" = model_matrix_result_prop %in% colnames(ptmp))
  x_regression_var_value = unlist(ptmp[,model_matrix_result_prop]) # model_matrix_result_prop is typically "condition"
  cols1 = x_regression_var_value == x_regression_var_value[1] # columns condition 1
  cols2 = !cols1                                              # columns condition 2
  # test for both NA and zeros
  # (compatible with eBayes/DEqMS/MSqRob ExpressionSets that have log2 intensities,
  #  and with MS-EmpiRe dataset with plain intensities sans log transformation that contain zeros for missing values)
  rows = rowSums(is.finite(x[,cols1,drop=FALSE] > 0) & x[,cols1,drop=FALSE] > 0) >= 2 &
    rowSums(is.finite(x[,cols2,drop=FALSE]) & x[,cols2,drop=FALSE] > 0) >= 2

  if(sum(rows) < 10) {
    append_log(sprintf("less than 10 %s have a non-zero value in both conditions: too little input data for DEA", ifelse(is_peptide, "peptides", "proteins")), type = "error")
  }

  # if there are any rows to remove, create a subset of the ExpressionSet + show warning (i.e. this should never occur during typical MS-DAP workflow that includes filtering)
  if(any(!rows)) {
    ftmp = Biobase::fData(eset)[rows,]
    y = x[rows,,drop=FALSE]
    eset = Biobase::ExpressionSet(
      assayData = y,
      featureData = Biobase::annotatedDataFrameFrom(y, byrow = T),
      protocolData = Biobase::annotatedDataFrameFrom(y, byrow = F)
    )
    Biobase::pData(eset) = ptmp
    Biobase::fData(eset) = ftmp

    append_log(sprintf("DEA requires at least 2 values per experimental condition; %d / %d %s were removed prior to applying %s\n(in typical MS-DAP workflows one would apply peptide-filtering upstream, e.g. with parameters in the analysis_quickstart() function)",
                       sum(!rows), length(rows), ifelse(is_peptide, "peptides", "proteins"), method_name), type = "warning")
  }

  return(eset)
}



#' Apply limma::eBayes() function to a protein-level ExpressionSet
#'
#' ref; PMID:25605792
#' ref; https://bioconductor.org/packages/release/bioc/html/limma.html
#'
#' @examples \dontrun{
#'   # Minimal example:
#'   # Let `m` be a log2 protein intensity matrix matrix where
#'   # columns are samples and rows proteins.
#'   # Let `eset` be a Biobase::ExpressionSet that contains `m`.
#'   # Let `m_groups` be a character vector that describes the
#'   # sample group/condition for each column in `m`.
#'   # e.g. m_groups = c("WT", "WT", "WT", "KO", "KO", "KO")
#'   model_matrix = stats::model.matrix(~m_groups)
#'   de_ebayes(eset, model_matrix, colnames(model_matrix)[2])
#' }
#' @param eset protein-level log2 intensity matrix stored as a Biobase ExpressionSet
#' @param model_matrix a `stats::model.matrix()` result that is supplied to `limma::lmFit()`
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
#' @export
de_ebayes = function(eset, model_matrix, model_matrix_result_prop) {
  start_time = Sys.time()

  # will throw an error if input data is invalid
  validate_eset_fitdata(eset, model_matrix, model_matrix_result_prop)
  # ensure there are at least 2 valid data points per condition / group @ model_matrix_result_prop (subset the expressionset otherwise)
  eset = subset_protein_eset_for_dea(eset, model_matrix_result_prop, is_peptide = FALSE, method_name = "eBayes")

  x = Biobase::exprs(eset)
  fit = suppressWarnings(limma::eBayes(limma::lmFit(x, model_matrix)))
  # !! sort.by="none" keeps the output table aligned with input matrix
  result = suppressMessages(limma::topTable(fit, number = nrow(x), coef = model_matrix_result_prop, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  # note; fit$coefficients and fit$stdev.unscaled matrices contain the intercept, don't use in the ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  result$effectsize = fit$coefficients[,model_matrix_result_prop] / sqrt(fit$s2.post)
  # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  result$standarderror = sqrt(fit$s2.post) * fit$stdev.unscaled[,model_matrix_result_prop]
  result$standarddeviation = sqrt(fit$s2.post)

  # convert from data.frame to a tibble that follows the column names/format we expect downstream
  result = as_tibble(result) %>%
    mutate(protein_id = rownames(x)) %>%
    select(protein_id, pvalue = P.Value, qvalue = adj.P.Val, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation, standarderror) %>%
    add_column(dea_algorithm = "ebayes")

  append_log_timestamp("eBayes", start_time)
  return(result)
}



#' Apply DEqMS to a protein-level ExpressionSet
#'
#' ref; PMID:32205417
#' ref; https://github.com/yafeng/DEqMS
#'
#' implementation follows example code from the vignette; https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html
#'
#' @param eset protein-level log2 intensity matrix stored as a Biobase ExpressionSet. Note that it must contain protein metadata column 'npep' that holds integer peptide counts
#' @param model_matrix a `stats::model.matrix()` result that is supplied to `limma::lmFit()`
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
#' @param doplot create a QC plot?
#' @export
de_deqms = function(eset, model_matrix, model_matrix_result_prop, doplot = FALSE) {
  start_time = Sys.time()

  # will throw an error if input data is invalid
  validate_eset_fitdata(eset, model_matrix, model_matrix_result_prop)
  # ensure there are at least 2 valid data points per condition / group @ model_matrix_result_prop (subset the expressionset otherwise)
  eset = subset_protein_eset_for_dea(eset, model_matrix_result_prop, is_peptide = FALSE, method_name = "DEqMS")

  # additional input validation; DEqMS requires peptide counts
  tmp = Biobase::fData(eset)$npep
  if(length(tmp) == 0 || !all(is.finite(tmp) & is.integer(tmp) & tmp > 0)) {
    append_log("parameter eset must contain a column 'npep' (i.e. Biobase::fData(eset)) with positive integer values", type = "error")
  }
  rm(tmp)

  # eBayes fit
  fit = suppressWarnings(limma::eBayes(limma::lmFit(Biobase::exprs(eset), model_matrix)))

  ### bugfix for DEqMS::spectraCounteBayes()
  # DEqMS will Loess fit log peptide counts versus log sigma^2
  # However, in some datasets/normalizations a subset of `fit$sigma` values from limma::eBayes() fit
  # Log transforming sigma^0 where sigma is zero results in a non-finite value, which crashes the loess fit
  # fix: replace zero sigma's with lowest non-zero sigma  (easy fix without having to update DEqMS code)
  # ? perhaps it's better if DEqMS would fit against eBayes' posterior estimates of sigma instead ?
  fit$sigma[!is.finite(fit$sigma) | fit$sigma <= 0] = min(fit$sigma[is.finite(fit$sigma) & fit$sigma > 0])

  ## peptide-per-protein counts
  # note that the ebayes fit results are aligned with the intensity-value-matrix we supplied upstream,
  # which in turn is aligned with the metadata in the ExpressionSet object. So we can just use its npep property as-is
  fit$count = Biobase::fData(eset)$npep

  # apply DEqMS, then overwrite the fit object with DEqMS's
  fit = suppressWarnings(DEqMS::spectraCounteBayes(fit))
  if(doplot) {
    suppressWarnings(DEqMS::VarianceBoxplot(fit, n=20, main = "DEqMS QC plot", xlab="#unique peptides per protein"))
  }

  ## analogous to our de_ebayes() implementation, extract results from the fit object but grab the DEqMS specific output columns where available (reference; DEqMS::outputResult() )
  # !! sort.by="none" keeps the output table aligned with input matrix
  result = suppressMessages(limma::topTable(fit, number = length(fit$count), coef = model_matrix_result_prop, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  # indexing by name might fail in the deqms output table for some model matrices !
  if(model_matrix_result_prop %in% colnames(fit$sca.p)) {
    result$pvalue = fit$sca.p[,model_matrix_result_prop]
    result$tstatistic = fit$sca.t[,model_matrix_result_prop]
  } else {
    result$pvalue = fit$sca.p[,match(model_matrix_result_prop, colnames(fit$coefficients))] # fit$sca.p might be unnamed -> use "standard" property of similar dimensions to find index
    result$tstatistic = fit$sca.t[,match(model_matrix_result_prop, colnames(fit$coefficients))]
  }
  result$qvalue = p.adjust(result$pvalue, method = "fdr")
  ## fit$coefficients and fit$stdev.unscaled contain the intercept, remove from ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  result$effectsize = fit$coefficients[,model_matrix_result_prop] / sqrt(fit$sca.postvar)
  # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  result$standarderror = sqrt(fit$sca.postvar) * fit$stdev.unscaled[,model_matrix_result_prop]
  result$standarddeviation = sqrt(fit$sca.postvar)

  ## create a result tibble that contains all columns required for downstream compatability with this pipeline; protein_id, pvalue, qvalue, foldchange.log2, dea_algorithm
  result = as_tibble(result) %>%
    mutate(protein_id = rownames(fit)) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic, standarddeviation, standarderror) %>%
    add_column(dea_algorithm = "deqms")

  append_log_timestamp("DEqMS", start_time)
  return(result)
}



#' MS-EmpiRe implementation, a wrapper function for msEmpiRe::de.ana()
#'
#' ref; PMID:31235637
#' ref; https://github.com/zimmerlab/MS-EmpiRe
#'
#' note; MS-EmpiRe crashes when all proteins have the same number of peptides (eg; filtering topn=1 or topn=3 and minpep=3)
#'
#' TODO: proper standarderror computation
#'
#' @param eset a Biobase ExpressionSet that contains the peptide log2 intensities
#' @param model_matrix a `stats::model.matrix()` result that is supplied to `limma::lmFit()`
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
#' @export
de_msempire = function(eset, model_matrix, model_matrix_result_prop) {
  start_time = Sys.time()

  # issue with MS-EmpiRe we previously reported, including code to reproduce and resolve this; https://github.com/zimmerlab/MS-EmpiRe/issues/10
  # bugfix: fix the random seed to ensure MS-EmpiRe results are always the same (given the exact same input data)
  set.seed(123)

  # will throw an error if input data is invalid
  validate_eset_fitdata(eset, model_matrix, model_matrix_result_prop)
  # MS-EmpiRe will crash if there are input values with all zeros
  # ensure there are at least 2 valid data points per condition / group @ model_matrix_result_prop (subset the expressionset otherwise)
  eset = subset_protein_eset_for_dea(eset, model_matrix_result_prop, is_peptide = FALSE, method_name = "MS-EmpiRe")

  # importantly, MS-EmpiRe expects the input data matrix to be NOT be log2 transformed
  x = Biobase::exprs(eset)
  x = 2^x
  x[!is.finite(x)] = 0
  Biobase::exprs(eset) = x

  # protein identifiers are expected in column "prot.id"
  Biobase::fData(eset)["prot.id"] = Biobase::fData(eset)[, "protein_id"]

  # experimental condition is expected in column "condition"
  # we here use the values from the design matrix @ desired output column
  Biobase::pData(eset)["condition"] = unname(model_matrix[,model_matrix_result_prop])

  # call msempire main function while silencing all of its output
  tmp = capture.output(result <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(eset))))
  append_log_timestamp("MS-EmpiRe", start_time)

  tibble(
    protein_id = as.character(result$prot.id), # enforce type to be robust against upstream R package changes (e.g. prot.id as a factor)
    pvalue = as.numeric(result$p.val),
    qvalue = as.numeric(result$p.adj),
    foldchange.log2 = as.numeric(result$log2FC),
    # update;
    # In MS-EmpiRe, 'prot.sd' is the standard deviation of 'prot.s' (which underwent shrinkage),
    # not of the protein foldchange. So here we compute the effectsize based on prot.s/prot.sd analogous to
    # MS-EmpiRe's protein p-value computation.
    #
    # However, there seems to be a bug in MS-EmpiRe where the sign of prot.s does not always agree with log2FC !
    # So this simple test fails on several datasets;
    # stopifnot( (result$log2FC >= 0) == (result$prot.s / result$prot.sd >= 0) )
    #
    # As a workaround, we compute absolute effectsize and then apply the sign of log2FC
    effectsize = abs(as.numeric(result$prot.s) / as.numeric(result$prot.sd)) * ifelse(foldchange.log2 < 0, -1, 1),
    #previous# effectsize = as.numeric(result$log2FC) / as.numeric(result$prot.sd),
    tstatistic = NA,
    # approximate standard deviation for the log2FC from the standardized effectsize
    standarddeviation = abs(foldchange.log2 / effectsize),
    #previous# standarddeviation = as.numeric(result$prot.sd),
    standarderror = NA,
    dea_algorithm = "msempire"
  ) %>%
    # if either log2FC or p.val is invalid, remove the protein from results
    filter(is.finite(foldchange.log2) & is.finite(pvalue))
}



#' Wrapper function for msqrob, using the implementation from the msqrobsum package
#'
#' note that this is the final version of this algorithm provided by its authors, and is no longer maintained; https://github.com/statOmics/msqrob
#' only difference with the GitHub version is that code included here has minor adaptions to situationally improve multithreading (i.e. speeds up computation a bit, on some systems)
#'
#' ref; PMID:26566788 PMID:32321741
#' ref; https://github.com/statOmics/msqrob
#'
#' @param eset a Biobase ExpressionSet that contains the peptide log2 intensities
#' @param model_matrix a `stats::model.matrix()` result that is supplied to `limma::lmFit()`
#' @param model_matrix_result_prop the column name in `model_matrix` that should be returned as the resulting coefficient estimated by `eBayes()`. In MS-DAP workflows this is typically "condition"
#' @param use_peptide_model if `TRUE`, apply msqrob. if `FALSE`, apply msqrobsum
#' @param random_variables optionally, an array of column names in the sample metadata table that should be used as additional regression variables
#' @export
de_msqrobsum_msqrob = function(eset, model_matrix, model_matrix_result_prop, use_peptide_model, random_variables = NULL) {
  start_time = Sys.time()

  # will throw an error if input data is invalid
  validate_eset_fitdata(eset, model_matrix, model_matrix_result_prop)
  if(length(use_peptide_model) != 1 || !use_peptide_model %in% c(TRUE, FALSE)) {
    append_log("use_peptide_model parameter must be either TRUE or FALSE", type = "error")
  }

  algorithm_name = ifelse(use_peptide_model, "msqrob", "msqrobsum")

  # ensure there are at least 2 valid data points per condition / group @ model_matrix_result_prop (subset the expressionset otherwise)
  eset = subset_protein_eset_for_dea(eset, model_matrix_result_prop, is_peptide = TRUE, method_name = algorithm_name)

  # just to be sure, fix random seed
  set.seed(123)


  #### below text quoted from MSqRobsum manual   @   https://github.com/statOmics/MSqRobSum/blob/master/vignettes/msqrobsum.Rmd
  # We can also use the `msqrobsum()` function to perform a MSqRob analysis on peptide intensities without first summarizing to protein summaries.
  # Because a protein can have intensities from multiple peptides and the intensities belonging to 1 peptide are correlated with eachother whe have to account for this in our model.
  # Previously we only had 1 protein summary per sample but now we have multiple peptide intensities per sample and these are also correlated with eachother.
  # Hence our model: `expression ~ (1|condition) + (1|sample) + (1|feature)`.
  # However some proteins will only have intensities from 1 peptide and the model fitting wil fail if we try to use the model above. For these proteins we should use te reduced model `expression ~ (1|condition)`. `msqrobsum()`
  # The  `formulas ` parameter accepts a vector of formulas. Model fitting with the first model will be attempted first but if that fails it tries the second model and so on.

  msnset = MSnbase::as.MSnSet.ExpressionSet(eset)
  if (use_peptide_model) {
    form_string = c(sprintf("expression ~ (1 | %s) + (1 | sample_id) + (1 | peptide_id)", model_matrix_result_prop), sprintf("expression ~ (1 | %s)", model_matrix_result_prop))
  } else {
    form_string = sprintf("expression ~ (1 | %s)", model_matrix_result_prop)
  }

  # strip forbidden random variables
  random_variables = setdiff(random_variables, c("peptide_id", "protein_id", "sample_id", model_matrix_result_prop))

  # compose updated formula if the user provided additional random variables
  if(length(random_variables) > 0) {
    form_string = sapply(form_string, function(x) paste(c(x, sprintf("(1 | %s)", random_variables)), collapse=" + "), USE.NAMES = FALSE)
  }

  # always add y~(1|condition) without any optional random_variables as a last-resort model if all other model specifications fail (due to lack of data) -->> then convert from string to formula/expression
  form_string = unique(c(form_string, sprintf("expression ~ (1 | %s)", model_matrix_result_prop)))
  form_eval = sapply(form_string, function(x) eval(parse(text=x)), USE.NAMES = FALSE)
  append_log(sprintf("%s linear regression formula%s; %s", algorithm_name, ifelse(length(form_string)>1, "s (these are prioritized. eg; if a model fit fails due to lack of data, the next formula is used)", ""), paste(form_string, collapse = "  ,  ")), type = "info")


  # MSqRob model as re-implemented in msqrobsum package, by original authors, and updated by us (only difference is computational performance)
  if (use_peptide_model) {
    result = suppressWarnings(msqrobsum(data = msnset, formulas = form_eval, group_vars = "protein_id", contrasts = model_matrix_result_prop, mode = "msqrob"))
    # result = suppressWarnings(msqrobsum(data = msnset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrob"))
    # our version 1, without user-specified random variables; form = c(expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id), expression ~ (1 | condition))
  } else {
    ### first use MSqRobSum to rollup to protein level, 'robust' approach appears to be an improvement over the traditional 'sum', then apply msqrob
    # line 273  @  https://github.com/statOmics/MSqRobSum/blob/97e22fddb9d6f1d3c29aafbae28c382148b9471d/vignettes/msqrobsum.Rmd#L273
    protset = suppressWarnings(suppressMessages(MSnbase::combineFeatures(msnset, fun = "robust", groupBy = Biobase::fData(msnset)$protein_id)))
    # line 313  @  https://github.com/statOmics/MSqRobSum/blob/97e22fddb9d6f1d3c29aafbae28c382148b9471d/vignettes/msqrobsum.Rmd#L313
    result = suppressWarnings(msqrobsum(data = protset, formulas = form_eval, group_vars = "protein_id", contrasts = model_matrix_result_prop, mode = "msqrobsum"))
    # result = suppressWarnings(msqrobsum(data = protset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrobsum"))
    ## our version 1, without user-specified random variables;
    #result = suppressWarnings(msqrobsum(data = protset, formulas = expression ~ (1 | condition), group_vars = "protein_id", contrasts = "condition", mode = "msqrobsum"))
  }

  result_unpacked = as_tibble(result) %>%
    dplyr::select(protein_id, contrasts) %>%
    tidyr::unnest(cols = contrasts)

  result_unpacked$sigma_post = result_unpacked$se / result_unpacked$sigma_contrast
  result_unpacked$effectsize = result_unpacked$logFC / result_unpacked$sigma_post

  # prepare output table in our standard format
  result = as_tibble(result_unpacked) %>%
    mutate(protein_id = as.character(protein_id),
           dea_algorithm = algorithm_name) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation = sigma_post, standarderror = se, dea_algorithm)

  append_log_timestamp(algorithm_name, start_time)
  return(result)
}



#' Perform linear regression with limma
#'
#' @description
#' Please refer to the vignette on custom limma models at the MS-DAP GitHub repository for a walk-through with examples.
#'
#' @param x result from `get_protein_matrix(... , include_npep = TRUE)`. Alternatively, a protein log2 intensity matrix (in this case, also provide the `npep` parameter if you want to apply DEqMS)
#' @param npep optionally, the number of peptides-per-protein for each row in input protein matrix. Only needed when providing a matrix for parameter `x` (i.e. you can leave this NULL when `x` is output from `get_protein_matrix(... , include_npep = TRUE)`)
#' @param model_matrix result from `stats::model.matrix()`
#' @param contrasts contrasts defined with `limma::makeContrasts()` that should be extracted from the linear regression model
#' @param limma_block_variable option passed to `limma::lmFit()` to fit block design. To skip (default), set NULL. Otherwise, this should be a character vector of the same length as nrow(model_matrix) so it describes a block/group per respective sample in the model/design matrix
#' @param ebayes_trend option passed to `limma::eBayes()`. Defaults to `FALSE` (same as limma default)
#' @param ebayes_robust option passed to `limma::eBayes()`. Defaults to `FALSE` (same as limma default)
#' @param deqms apply DEqMS? Defaults to `TRUE`. If `FALSE`, returns `limma::eBayes()` result as-is
#' @param return_table if `TRUE`, the default, returns a table with protein_id, log2fc, pvalue, adjusted pvalue, contrast (i.e. this is the output from MS-DAP function `limma_fit_extract_stats()`). If `FALSE`, returns the limma fit object as-is
#' @export
limma_wrapper = function(
    x = protein_data,
    npep = NULL,
    model_matrix,
    contrasts,
    limma_block_variable = NULL,
    ebayes_trend = FALSE,
    ebayes_robust = FALSE,
    deqms = TRUE,
    return_table = TRUE
) {

  mat = NULL

  ### input validation
  # test that we got the proper data formats (matrix and vector)
  if(!is.matrix(x)) {
    if(!(is.list(x) && all(c("matrix", "npep") %in% names(x)) && is.matrix(x$matrix) &&
         (is.vector(x$npep) || is.array(x$npep)) && is.numeric(x$npep) && length(x$npep) == nrow(x$matrix))) {
      append_log("x parameter should either be a matrix, or a list that contains a 'matrix' and 'npep' parameter (where npep is an integer vector of same length as nrow(x$matrix)). Typically this is the output from MS-DAP function: get_protein_matrix(... , include_npep = TRUE)", type = "error")
    }
    if(!is.null(npep)) {
      append_log("when parameter x is a list (i.e. input from function get_protein_matrix(... , include_npep = TRUE)), leave the npep parameter empty", type = "error")
    }

    mat = x$matrix
    npep = as.vector(x$npep)
  } else {
    if(!((is.vector(npep) || is.array(npep)) && is.numeric(npep) && length(npep) == nrow(x))) {
      append_log("when parameter x is a matrix, the npep parameter should be an integer vector of same length as nrow(x)", type = "error")
    }
    mat = x
    npep = as.vector(npep)
  }

  # test that the matrix is numeric and that npep only contains finite values
  if(mode(mat) != "numeric") {
    append_log("matrix must be numeric (you can check with the mode() function)", type = "error")
  }
  if(!all(is.finite(npep) & npep > 0)) {
    append_log("npep must only contain finite integers, larger than zero", type = "error")
  }

  # finally, check simple parameters
  if(!is.matrix(model_matrix)) {
    append_log("parameter model_matrix must be a matrix, i.e. output from stats::model.matrix()", type = "error")
  }
  if(!is.matrix(contrasts)) {
    append_log("parameter contrasts must be a matrix, i.e. output from limma::makeContrasts()", type = "error")
  }
  if(!is.null(limma_block_variable) && (
    anyNA(limma_block_variable) || !is.character(limma_block_variable) ||
    any(limma_block_variable == "") || length(limma_block_variable) != nrow(model_matrix))
  ) {
    append_log("parameter limma_block_variable must be a character vector of the same length as nrow(model_matrix) so it describes a block/group per respective sample in the model/design matrix. Empty strings and NA not allowed", type = "error")
  }
  if(length(ebayes_trend) != 1 || !ebayes_trend %in% c(TRUE, FALSE)) {
    append_log("parameter ebayes_trend must be either TRUE or FALSE", type = "error")
  }
  if(length(ebayes_robust) != 1 || !ebayes_robust %in% c(TRUE, FALSE)) {
    append_log("parameter ebayes_robust must be either TRUE or FALSE", type = "error")
  }
  if(length(deqms) != 1 || !deqms %in% c(TRUE, FALSE)) {
    append_log("parameter deqms must be either TRUE or FALSE", type = "error")
  }
  if(length(return_table) != 1 || !return_table %in% c(TRUE, FALSE)) {
    append_log("parameter return_table must be either TRUE or FALSE", type = "error")
  }
  ### input validation


  ### prepare data
  # enforce integer type. use ceiling so 0.1 becomes 1 (i.e. we never have zeros)
  npep = as.integer(ceiling(npep))


  ### apply limma
  # adapted from limma documentation; https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf , section 9.7 "Multi-level Experiments"
  fit = NULL
  if(is.null(limma_block_variable)) {
    fit = limma::lmFit(mat, model_matrix)
  } else {
    corfit = limma::duplicateCorrelation(mat, model_matrix, block = limma_block_variable)
    append_log(sprintf('applying a limma block design: %.3f correlation between samples within the same "block"\n', corfit$consensus), type = "info")
    fit = limma::lmFit(mat, model_matrix, block = limma_block_variable, correlation = corfit$consensus)
  }

  fit2 = limma::contrasts.fit(fit, contrasts)
  fit2 = limma::eBayes(fit2, trend = ebayes_trend, robust = ebayes_robust)


  ### apply DEqMS, analogous to msdap::de_deqms()
  if(deqms) {
    # bugfix for DEqMS::spectraCounteBayes(), see
    fit2$sigma[!is.finite(fit2$sigma) | fit2$sigma <= 0] = min(fit2$sigma[is.finite(fit2$sigma) & fit2$sigma > 0])
    # add peptide-per-protein counts to the limma fit object, as required for DEqMS
    fit2$count = npep
    fit2 = suppressWarnings(DEqMS::spectraCounteBayes(fit2))
    append_log("applying limma eBayes(), followed by DEqMS (to adjust protein confidences according to their respective number of peptides)", type = "info")
  } else {
    append_log("returning limma eBayes() results as-is (i.e. no DEqMS post-hoc correction)", type = "info")
  }


  if(return_table) {
    return(limma_fit_extract_stats(fit2))
  } else {
    return(fit2)
  }
}



#' extract limma fit results
#'
#' @param fit result from `limma::lmFit()`
#' @param contrast_name name of a specific contrast(s) to extract. Default value, `NULL`, will return results from all contrasts
limma_fit_extract_stats = function(fit, contrast_name = NULL) {
  all_contrasts = setdiff(colnames(fit$coefficients), c("intercept", "Intercept", "(intercept)", "(Intercept"))

  # naive check for limma result
  if(!is.list(fit) || !all(c("coefficients", "sigma", "s2.post", "stdev.unscaled") %in% names(fit)) ) {
    append_log("fit parameter must be the result from limma::eBayes()  (e.g. obtained from MS-DAP function limma_wrapper() )", type = "error")
  }
  if(is.null(contrast_name)) {
    contrast_name = all_contrasts
  } else {
    if(anyNA(all_contrasts) || !is.character(contrast_name) || any(!contrast_name %in% all_contrasts)) {
      append_log(paste("some values in contrast_name are not found in the provided fit parameter. Available contrasts in fit:", paste(all_contrasts, collapse = ",")), type = "error")
    }
  }

  result = NULL
  for(contr in contrast_name) {
    x = tibble::as_tibble(limma::topTable(fit, number = Inf, coef = contr, adjust.method = "fdr", sort.by = "none", confint = TRUE))
    x$protein_id = rownames(fit)
    x$contrast = contr
    x$foldchange.log2 = x$logFC

    # test if the deqms R package was applied to this limma::eBayes() fit/result
    is_deqms = "sca.p" %in% names(fit)

    if(is_deqms) {
      ## DEqMS code analogous to msdap:::de_deqms()
      # indexing by name might fail in the deqms output table for some model matrices !
      if(contr %in% colnames(fit$sca.p)) {
        x$pvalue = fit$sca.p[,contr]
      } else {
        x$pvalue = fit$sca.p[,match(contr, colnames(fit$coefficients))] # fit$sca.p might be unnamed -> use "standard" property of similar dimensions to find index
      }
      x$qvalue = p.adjust(x$pvalue, method = "fdr") # in this loop, we subsetted a specific contrast so can safely call p.adjust()
      # analogous to standard ebayes fit @ below else block, but note that we use "sca" parameters added by DEqMS !
      x$effectsize = fit$coefficients[,contr] / sqrt(fit$sca.postvar)
      x$standarderror = sqrt(fit$sca.postvar) * fit$stdev.unscaled[,contr]
      x$standarddeviation = sqrt(fit$sca.postvar)
    } else {
      ## eBayes code analogous to msdap:::de_ebayes()
      x$pvalue = x$P.Value
      x$qvalue = x$adj.P.Val
      x$effectsize = fit$coefficients[,contr] / sqrt(fit$s2.post)
      x$standarderror = sqrt(fit$s2.post) * fit$stdev.unscaled[,contr]
      x$standarddeviation = sqrt(fit$s2.post)
    }

    result = bind_rows(
      result,
      x %>% select(protein_id, contrast, foldchange.log2, effectsize, standarddeviation, standarderror, pvalue, qvalue) %>% mutate_all(unname)
    )
  }

  return(result)
}
