
#' apply limma::eBayes() function to a protein-level ExpressionSet
#'
#' ref; PMID:25605792
#' ref; https://bioconductor.org/packages/release/bioc/html/limma.html
#'
#' @param eset_proteins protein-level dataset stored as a Biobase ExpressionSet
#' @param input_intensities_are_log2 boolean indicating whether the ExpressionSet's intensity values are already log2 scaled (default: TRUE)
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms in each statistical contrast tested downstream
#'
#' @importFrom Biobase exprs pData
#' @importFrom limma eBayes topTable lmFit
#' @export
de_ebayes = function(eset_proteins, input_intensities_are_log2 = TRUE, random_variables = NULL) {
  start_time = Sys.time()
  if(length(random_variables) > 0 && (!is.vector(random_variables) || any(is.na(random_variables) | !is.character(random_variables))) ) {
    append_log(paste("random_variables must either be NULL or character vector", paste(random_variables, collapse=", ")), type = "error")
  }

  # transform to log2 if input data is non-log
  if (!input_intensities_are_log2) {
    x = log2(Biobase::exprs(eset_proteins))
    x[!is.finite(x)] = NA
    Biobase::exprs(eset_proteins) = x
  }

  # data.frame with sample metadata for those random variables that are to be tested in limma::eBayes
  random_variables = unique(c("condition", random_variables)) # use unique just in case the user erronously passed the condition as 'additional random variable'
  df_metadata = Biobase::pData(eset_proteins)

  if(!all(random_variables %in% colnames(df_metadata))) {
    append_log(paste("eset_proteins ExpressionSet pData() must contain columns 'condition' and all additional random variables requested in `random_variables` parameter. missing:", paste(setdiff(random_variables, colnames(df_metadata)), collapse=", ")), type = "error")
  }

  # input validation, then apply limma::ebayes()
  x = Biobase::exprs(eset_proteins)
  df_metadata = df_metadata[ , random_variables, drop=F]
  fit = de_ebayes_fit(x, random_variables = df_metadata)

  # !! sort.by="none" keeps the output table aligned with input matrix
  # always extract the first regression variable (intercept has index 1, so the column number we need is 2)
  result = suppressMessages(limma::topTable(fit, number = nrow(x), coef = 2, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  # note; fit$coefficients and fit$stdev.unscaled matrices contain the intercept, don't use in the ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  result$effectsize = fit$coefficients[,2] / sqrt(fit$s2.post)
  # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  result$standarderror = sqrt(fit$s2.post) * fit$stdev.unscaled[,2]
  result$standarddeviation = sqrt(fit$s2.post)

  # convert from data.frame to a tibble that follows the column names/format we expect downstream
  result = as_tibble(result) %>%
    mutate(protein_id = rownames(x)) %>%
    select(protein_id, pvalue = P.Value, qvalue = adj.P.Val, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation, standarderror) %>%
    add_column(dea_algorithm = "ebayes")

  # note; we added condition to random variables above, strip it from result table for consistency within this R package ('condition' is always assumed to be a regression variable)
  if(length(random_variables) > 1) { # ! test larger than 1
    result = result %>% add_column(contrast_ranvars = paste(setdiff(random_variables, "condition"), collapse = ", "))
  }
  append_log_timestamp("eBayes", start_time)
  return(result)
}



#' apply DEqMS to a protein-level ExpressionSet
#'
#' ref; PMID:32205417
#' ref; https://github.com/yafeng/DEqMS
#'
#' implementation follows example code from the vignette; https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html
#'
#' @param eset_proteins protein-level dataset stored as a Biobase ExpressionSet. Note that it must contain protein metadata column 'npep' that holds integer peptide counts and sample metadata column 'condition'
#' @param input_intensities_are_log2 boolean indicating whether the ExpressionSet's intensity values are already log2 scaled (default: TRUE)
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms in each statistical contrast tested downstream
#' @param doplot create a QC plot?
#'
#' @importFrom Biobase exprs pData
#' @importFrom limma eBayes topTable lmFit
#' @importFrom DEqMS spectraCounteBayes
#' @export
de_deqms = function(eset_proteins, input_intensities_are_log2 = TRUE, random_variables = NULL, doplot = FALSE) {
  start_time = Sys.time()

  ### input validation

  if(length(random_variables) > 0 && (!is.vector(random_variables) || any(is.na(random_variables) | !is.character(random_variables))) ) {
    append_log(paste("random_variables must either be NULL or character vector", paste(random_variables, collapse=", ")), type = "error")
  }

  if(length(Biobase::pData(eset_proteins)) == 0 || !is.data.frame(Biobase::pData(eset_proteins))) {
    append_log("eset_proteins ExpressionSet pData must be a data.frame", type = "error")
  }

  if(length(Biobase::fData(eset_proteins)) == 0 || !is.data.frame(Biobase::fData(eset_proteins))) {
    append_log("eset_proteins ExpressionSet fData must be a data.frame", type = "error")
  }

  if(length(Biobase::pData(eset_proteins)$condition) == 0) {
    append_log("eset_proteins ExpressionSet pData must contain column 'condition'", type = "error")
  }

  tmp = Biobase::fData(eset_proteins)$npep
  if(length(tmp) == 0 || !all(is.finite(tmp) & is.integer(tmp) & tmp > 0)) {
    append_log("eset_proteins ExpressionSet fData must contain column 'npep' with positive integer values", type = "error")
  }
  rm(tmp)

  # data.frame with sample metadata for those random variables that are to be tested in limma::eBayes
  random_variables = unique(c("condition", random_variables)) # use unique just in case the user erronously passed the condition as 'additional random variable'
  df_metadata = Biobase::pData(eset_proteins)

  if(!all(random_variables %in% colnames(df_metadata))) {
    append_log(paste("eset_proteins ExpressionSet pData() must contain columns 'condition' and all additional random variables requested in `random_variables` parameter. missing:", paste(setdiff(random_variables, colnames(df_metadata)), collapse=", ")), type = "error")
  }

  ### input validation done


  # transform to log2 if input data is non-log
  if(!input_intensities_are_log2) {
    x = log2(Biobase::exprs(eset_proteins))
    x[!is.finite(x)] = NA
    Biobase::exprs(eset_proteins) = x
    rm(x)
    gc()
  }

  # input validation, then apply limma::ebayes()
  df_metadata = df_metadata[ , random_variables, drop=F]
  eset_proteins__protein_id = rownames(Biobase::exprs(eset_proteins))
  fit = de_ebayes_fit(Biobase::exprs(eset_proteins), random_variables = df_metadata)


  #### DEqMS unique code, everything else in this function is analogous to de_ebayes implementation ####

  ### bugfix for DEqMS::spectraCounteBayes()
  # DEqMS will Loess fit log peptide counts versus log sigma^2
  # However, in some datasets/normalizations a subset of `fit$sigma` values from limma::eBayes() fit
  # Log transforming sigma^0 where sigma is zero results in a non-finite value, which crashes the loess fit
  # fix: replace zero sigma's with lowest non-zero sigma  (easy fix without having to update DEqMS code)
  # ? perhaps it's better if DEqMS would fit against eBayes' posterior estimates of sigma instead ?
  fit$sigma[fit$sigma <= 0] = min(fit$sigma[fit$sigma > 0])

  ## peptide-per-protein counts
  # note that the ebayes fit results are aligned with the intensity-value-matrix we supplied upstream,
  # which in turn is aligned with the metadata in the ExpressionSet object. So we can just use its npep property as-is
  fit$count = Biobase::fData(eset_proteins)$npep

  # apply DEqMS, then overwrite the fit object with DEqMS's
  fit = suppressWarnings(DEqMS::spectraCounteBayes(fit))
  if(doplot) {
    suppressWarnings(DEqMS::VarianceBoxplot(fit, n=20, main = "DEqMS QC plot", xlab="#unique peptides per protein"))
  }

  ## analogous to our de_ebayes() implementation, extract results from the fit object but grab the DEqMS specific output columns where available (reference; DEqMS::outputResult() )
  coef_col = 2
  # !! sort.by="none" keeps the output table aligned with input matrix
  result = suppressMessages(limma::topTable(fit, number = length(eset_proteins__protein_id), coef = coef_col, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  stopifnot(rownames(result) == rownames(fit)) # because we didn't sort, tables are aligned
  result$tstatistic = fit$sca.t[,coef_col]
  result$pvalue = fit$sca.p[,coef_col]
  result$qvalue = p.adjust(result$pvalue, method = "fdr")
  ## fit$coefficients and fit$stdev.unscaled contain the intercept, remove from ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  result$effectsize = fit$coefficients[,coef_col] / sqrt(fit$sca.postvar)
  # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  result$standarderror = sqrt(fit$sca.postvar) * fit$stdev.unscaled[,coef_col]
  result$standarddeviation = sqrt(fit$sca.postvar)

  ## create a result tibble that contains all columns required for downstream compatability with this pipeline; protein_id, pvalue, qvalue, foldchange.log2, dea_algorithm
  result = as_tibble(result) %>%
    mutate(protein_id = eset_proteins__protein_id) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic, standarddeviation, standarderror) %>%
    add_column(dea_algorithm = "deqms")
  #### DEqMS unique code, everything else in this function is analogous to de_ebayes implementation ####


  # note; we added condition to random variables above, strip it from result table for consistency within this R package ('condition' is always assumed to be a regression variable)
  if(length(random_variables) > 1) { # ! test larger than 1
    result = result %>% add_column(contrast_ranvars = paste(setdiff(random_variables, "condition"), collapse = ", "))
  }
  append_log_timestamp("DEqMS", start_time)
  return(result)
}



#' Wrapper function for limma::ebayes()
#'
#' rownames of the data matrix (log2 intensities) must be the protein_id for downstream compatibility
#'
#' @param x log2 transformed protein intensity matrix
#' @param random_variables data.frame where each column describes a random variable (and each row matches a column in x) or a vector assumed to describe 1 random variable, sample group/condition for each column in x. Importantly, statistical results are returned only for the first random variable (so in most use-cases the sample group/condition should be described in the first column)
#'
#' @importFrom stats model.matrix
#' @importFrom limma eBayes lmFit
de_ebayes_fit = function(x, random_variables) {
  ## validate input; x must be a numeric matrix with rownames
  if(!is.matrix(x) || mode(x) != "numeric") {
    append_log("x must be a numerical matrix", type = "error")
  }
  if(length(rownames(x)) == 0) {
    append_log("input matrix must have rownames (these are assumed to be protein identifiers)", type = "error")
  }

  ## validate input; variables must not contain any NA and align with matrix x
  if(any(is.na(random_variables))) {
    append_log("random_variables must not contain any NA values", type = "error")
  }
  if(is.matrix(random_variables)) {
    random_variables = data.frame(random_variables, stringsAsFactors = F)
  }

  # if the supplied matrix/data.frame of random variables has any rownames, these must align with colnames of x
  if(length(rownames(random_variables)) == ncol(x) && !(
    all( rownames(random_variables) == as.character(1:nrow(random_variables)) ) ||
    all( rownames(random_variables) == colnames(x) )
  )) {
    append_log("input matrix x and random_variables seem misaligned; the column names of x versus rownames of variables do not match", type = "error")
  }

  if(is.vector(random_variables) || is.array(random_variables)) {
    random_variables = data.frame(group=random_variables, stringsAsFactors = F)
  }
  if(!is.data.frame(random_variables)) {
    append_log("random_variables must be an array, matrix or data.frame", type = "error")
  }
  if(nrow(random_variables) != ncol(x)) {
    append_log("number of columns in x (samples) must match the metadata (random_variables) provided for regression", type = "error")
  }


  # ensure all random_variables are factors; mutate each column, convert to factor
  # simply calling as.factor won't do as it orders factor levels alphabetically, while we want to retain the sorting from input data
  random_variables = random_variables %>% mutate_all(function(x) factor(x, levels=unique(x)) )

  # design matrix; simply use all columns in the user-provided metadata table
  contr_design = stats::model.matrix(~ . , data = random_variables)
  ## version 1
  # group_by_cols = match(group_by_cols, unique(group_by_cols)) - 1
  # contr_design = stats::model.matrix(~group_by_cols)

  return(suppressWarnings(limma::eBayes(limma::lmFit(x, contr_design))))

  ### full implementation
  # fit = suppressWarnings(limma::eBayes(limma::lmFit(x, contr_design)))
  # # !! sort.by="none" keeps the output table aligned with input matrix
  # # always extract the first regression variable (intercept has index 1, so the column number we need is 2)
  # result = suppressMessages(limma::topTable(fit, number = nrow(x), coef = 2, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  # ## fit$coefficients and fit$stdev.unscaled contain the intercept, remove from ES and SE computation
  # # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  # # note;
  # result$effectsize = fit$coefficients[,2] / sqrt(fit$s2.post)
  # # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  # result$standarderror = sqrt(fit$s2.post) * fit$stdev.unscaled[,2]
  # result$standarddeviation = sqrt(fit$s2.post)
  #
  # append_log_timestamp("eBayes", start_time)
  # return(as_tibble(result) %>%
  #          mutate(protein_id = rownames(x)) %>%
  #          select(protein_id, pvalue = P.Value, qvalue = adj.P.Val, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation, standarderror))
}



#' MS-EmpiRe implementation, a wrapper function for msEmpiRe::de.ana()
#'
#' ref; PMID:31235637
#' ref; https://github.com/zimmerlab/MS-EmpiRe
#'
#' input data should NOT be log2 transformed (if it is, make sure to set parameter `input_intensities_are_log2` so we can revert this prior to passing data to MS-EmpiRe)
#'
#' note; MS-EmpiRe crashes when all proteins have the same number of peptides (eg; filtering topn=1 or topn=3 and minpep=3)
#'
#' TODO: proper standarderror computation
#'
#' @param eset a Biobase ExpressionSet that contains the peptide intensities. required attributes; fData(): protein_id and pData(): condition
#' @param input_intensities_are_log2 whether the provided intensities are already log2 transformed
#'
#' @importFrom Biobase exprs pData fData
#' @importFrom msEmpiRe de.ana
#' @export
de_msempire = function(eset, input_intensities_are_log2 = F) {
  if(!("protein_id" %in% colnames(Biobase::fData(eset)))) {
    append_log("ExpressionSet fData() must contain column 'protein_id'", type = "error")
  }
  if(!("condition" %in% colnames(Biobase::pData(eset)))) {
    append_log("ExpressionSet pData() must contain columns 'condition'", type = "error")
  }

  # issue with MS-EmpiRe we previously reported, including code to reproduce and resolve this; https://github.com/zimmerlab/MS-EmpiRe/issues/10
  # bugfix: fix the random seed to ensure MS-EmpiRe results are always the same (given the exact same input data)
  set.seed(123)

  start_time = Sys.time()
  # ms-empire expects non-logtransformed data without NA values
  x = Biobase::exprs(eset)
  if (input_intensities_are_log2) {
    x = 2^x
  }
  x[!is.finite(x)] = 0
  Biobase::exprs(eset) = x

  # protein identifiers are expected in this column
  Biobase::fData(eset)["prot.id"] = Biobase::fData(eset)[, "protein_id"]
  # call msempire main function while silencing all of its output
  capture.output(result <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(eset))))
  append_log_timestamp("MS-EmpiRe", start_time)

  return(
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
  )
}



#' Wrapper function for msqrob, using the implementation from the msqrobsum package
#'
#' note that this is the final version of this algorithm provided by its authors, and is no longer maintained; https://github.com/statOmics/msqrob
#' only difference with the GitHub version is that code included here has minor adaptions to situationally improve multithreading (i.e. speeds up computation a bit, on some systems)
#'
#' ref; PMID:26566788 PMID:32321741
#' ref; https://github.com/statOmics/msqrob
#'
#' @param eset a Biobase ExpressionSet that contains the peptide intensities. required attributes; fData(): protein_id and pData(): condition
#' @param eset_proteins only required if `log2fc_without_shrinkage == TRUE`
#' @param use_peptide_model if true, apply msqrob. if false, apply msqrobsum
#' @param input_intensities_are_log2 whether the provided intensities are already log2 transformed
#' @param protein_rollup_robust for msqrobsum analyses, whether to use 'robust' protein rollup (TRUE) or traditional peptide intensity summation (FALSE)
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms
#' @param log2fc_without_shrinkage boolean parameter, if FALSE (default) returns msqrob estimated log2 foldchanges. if TRUE, uses `protein_foldchange_from_ebayes()` instead (affecting log2fc and effectsize results)
#'
#' @importFrom Biobase exprs pData fData
#' @importFrom MSnbase as.MSnSet.ExpressionSet combineFeatures
#' @export
de_msqrobsum_msqrob = function(eset, eset_proteins = NULL, use_peptide_model = T, input_intensities_are_log2 = T, protein_rollup_robust = T, random_variables = NULL, log2fc_without_shrinkage = FALSE) {
  start_time = Sys.time()
  if(log2fc_without_shrinkage && is.null(eset_proteins)) stop("if log2fc_without_shrinkage is set, provide protein eset")

  # strip forbidden random variables and compose updated formula
  random_variables = setdiff(random_variables, c("peptide_id", "protein_id", "sample_id", "condition"))

  if(!("protein_id" %in% colnames(Biobase::fData(eset)))) {
    append_log("ExpressionSet fData() must contain column 'protein_id'", type = "error")
  }
  if(use_peptide_model && !("peptide_id" %in% colnames(Biobase::fData(eset)))) {
    append_log("ExpressionSet fData() must contain column 'peptide_id' (when using peptide-level msqrob model)", type = "error")
  }
  if(!all(c("condition", "sample_id") %in% colnames(Biobase::pData(eset)))) {
    append_log("ExpressionSet pData() must contain columns 'condition' and 'sample_id'", type = "error")
  }
  if(length(random_variables) > 0 && !all(random_variables %in% colnames(Biobase::pData(eset)))) {
    append_log("ExpressionSet pData() must contain columns matching all random_variables", type = "error")
  }

  # just to be sure, fix random seed
  set.seed(123)


  # transform to log2 if input data is non-log
  if (!input_intensities_are_log2) {
    x = log2(Biobase::exprs(eset))
    x[!is.finite(x)] = NA
    Biobase::exprs(eset) = x
  }


  #### below text quoted from MSqRobsum manual   @   https://github.com/statOmics/MSqRobSum/blob/master/vignettes/msqrobsum.Rmd
  # We can also use the `msqrobsum()` function to perform a MSqRob analysis on peptide intensities without first summarizing to protein summaries.
  # Because a protein can have intensities from multiple peptides and the intensities belonging to 1 peptide are correlated with eachother whe have to account for this in our model.
  # Previously we only had 1 protein summary per sample but now we have multiple peptide intensities per sample and these are also correlated with eachother.
  # Hence our model: `expression ~ (1|condition) + (1|sample) + (1|feature)`.
  # However some proteins will only have intensities from 1 peptide and the model fitting wil fail if we try to use the model above. For these proteins we should use te reduced model `expression ~ (1|condition)`. `msqrobsum()`
  # The  `formulas ` parameter accepts a vector of formulas. Model fitting with the first model will be attempted first but if that fails it tries the second model and so on.

  algorithm_name = "msqrob"
  msnset = MSnbase::as.MSnSet.ExpressionSet(eset)
  if (use_peptide_model) {
    form_string = c("expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id)", "expression ~ (1 | condition)")
    # form_string = "expression ~ (1 | peptide_id) + (1 | sample_id) + (1 | condition)"
  } else {
    form_string = "expression ~ (1 | condition)"
    algorithm_name = "msqrobsum"
  }

  # compose updated formula if the user provided additional random variables
  if(length(random_variables) > 0) {
    form_string = sapply(form_string, function(x) paste(c(x, sprintf("(1 | %s)", random_variables)), collapse=" + "), USE.NAMES = FALSE)
  }

  # always add y~(1|condition) without any optional random_variables as a last-resort model if all other model specifications fail (due to lack of data) -->> then convert from string to formula/expression
  form_string = unique(c(form_string, "expression ~ (1 | condition)"))
  form_eval = sapply(form_string, function(x) eval(parse(text=x)), USE.NAMES = FALSE)
  append_log(sprintf("%s linear regression formula%s; %s", algorithm_name, ifelse(length(form_string)>1, "s (these are prioritized. eg; if a model fit fails due to lack of data, the next formula is used)", ""), paste(form_string, collapse = "  ,  ")), type = "info")

  # MSqRob model as re-implemented in msqrobsum package, by original authors, and updated by us (only difference is computational performance)
  if (use_peptide_model) {
    result = suppressWarnings(msqrobsum(data = msnset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrob"))
    # result = suppressWarnings(msqrobsum(data = msnset, formulas = c(eval(parse(text=form_string)), expression ~ (1 | condition)), group_vars = "protein_id", contrasts = "condition", mode = "msqrob"))
    # our version 1, without user-specified random variables; form = c(expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id), expression ~ (1 | condition))
  } else {
    ### first use MSqRobSum to rollup to protein level, 'robust' approach appears to be an improvement over the traditional 'sum', then apply msqrob
    # line 273  @  https://github.com/statOmics/MSqRobSum/blob/97e22fddb9d6f1d3c29aafbae28c382148b9471d/vignettes/msqrobsum.Rmd#L273
    protset = suppressWarnings(suppressMessages(MSnbase::combineFeatures(msnset, fun = ifelse(protein_rollup_robust, "robust", "sum"), groupBy = Biobase::fData(msnset)$protein_id)))
    # line 313  @  https://github.com/statOmics/MSqRobSum/blob/97e22fddb9d6f1d3c29aafbae28c382148b9471d/vignettes/msqrobsum.Rmd#L313
    result = suppressWarnings(msqrobsum(data = protset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrobsum"))
    ## our version 1, without user-specified random variables;
    #result = suppressWarnings(msqrobsum(data = protset, formulas = expression ~ (1 | condition), group_vars = "protein_id", contrasts = "condition", mode = "msqrobsum"))
  }

  result_unpacked = as_tibble(result) %>%
    dplyr::select(protein_id, contrasts) %>%
    tidyr::unnest(cols = contrasts)

  # test: optional log2fc estimates without shrinkage
  if(log2fc_without_shrinkage) {
    df_fc = protein_foldchange_from_ebayes(eset_proteins, input_intensities_are_log2 = TRUE, random_variables = random_variables)
    i = match(as.character(result_unpacked$protein_id), df_fc$protein_id)
    stopifnot(!is.na(i))
    result_unpacked$logFC = df_fc$foldchange.log2[i]
    algorithm_name = paste0(algorithm_name, "_fc")
  }

  result_unpacked$sigma_post = result_unpacked$se / result_unpacked$sigma_contrast
  result_unpacked$effectsize = result_unpacked$logFC / result_unpacked$sigma_post

  # prepare output table in our standard format
  result = as_tibble(result_unpacked) %>%
    mutate(protein_id = as.character(protein_id),
           dea_algorithm = algorithm_name) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation = sigma_post, standarderror = se, dea_algorithm)

  if(length(random_variables) > 0) {
    result = result %>% add_column(contrast_ranvars = paste(random_variables, collapse = ", "))
  }
  append_log_timestamp(algorithm_name, start_time)
  return(result)
}



#' Compute protein log2 foldchanges from protein-level eBayes (i.e. simple, no shrinkage)
#'
#' implementation analogous to `de_ebayes()`
#'
#' @param eset a Biobase ExpressionSet that contains the protein intensities. required attributes; fData(): protein_id and pData(): condition
#' @param input_intensities_are_log2 whether the provided intensities are already log2 transformed
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms
#' @return data.frame with columns protein_id and foldchange.log2
protein_foldchange_from_ebayes = function(eset, input_intensities_are_log2 = TRUE, random_variables = NULL) {
  if(length(random_variables) > 0 && (!is.vector(random_variables) ||
                                       any(is.na(random_variables) | !is.character(random_variables)))) {
    append_log(paste("random_variables must either be NULL or character vector",
                     paste(random_variables, collapse = ", ")), type = "error")
  }
  if(!input_intensities_are_log2) {
    x = log2(Biobase::exprs(eset))
    x[!is.finite(x)] = NA
    Biobase::exprs(eset) = x
  }
  random_variables = unique(c("condition", random_variables))
  df_metadata = Biobase::pData(eset)
  if(!all(random_variables %in% colnames(df_metadata))) {
    append_log(paste("eset ExpressionSet pData() must contain columns 'condition' and all additional random variables requested in `random_variables` parameter. missing:",
                     paste(setdiff(random_variables, colnames(df_metadata)),
                           collapse = ", ")), type = "error")
  }
  x = Biobase::exprs(eset)
  df_metadata = df_metadata[, random_variables, drop = F]

  # re-use our limma wrapper to apply eBayes()
  fit = de_ebayes_fit(x, random_variables = df_metadata)

  ## analogous to our de_ebayes() implementation, extract results from the fit object
  coef_col = 2
  data.frame(protein_id = rownames(x), foldchange.log2 = fit$coefficients[,coef_col])
}
