
#' placeholder title
#' @param peptides todo
#' @param samples todo
#' @param eset_peptides todo
#' @param eset_proteins todo
#' @param input_intensities_are_log2 todo
#' @param random_variables todo
#'
#' @importFrom Biobase exprs pData
de_interface_ebayes = function(peptides=NULL, samples=NULL, eset_peptides=NULL, eset_proteins=NULL, input_intensities_are_log2 = T, random_variables = NULL) {
  ## one can either use the peptide and samples tibbles, or a pre-constructed peptide/protein-level ExpressionSet
  ## in this example code, we directly work with the protein ExpressionSet and apply the eBayes algorithm @ limma

  if(length(eset_proteins) == 0) {
    append_log("current implementation of this function only supports input data through 'eset_proteins' parameters; a protein-level ExpressionSet must be provided", type = "error")
  }
  if(!("condition" %in% colnames(Biobase::pData(eset_proteins)))) {
    append_log("eset_proteins ExpressionSet pData() must contain columns 'condition'", type = "error")
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

  df_metadata = df_metadata[ , random_variables, drop=F]
  result = de_ebayes(x = Biobase::exprs(eset_proteins), random_variables = df_metadata) %>% add_column(algo_de = "ebayes")

  # note; we added condition to random variables above, strip it from result table for consistency within this R package ('condition' is always assumed to be a regression)
  if(length(random_variables) > 1) { # ! test larger than 1
    result = result %>% add_column(contrast_ranvars = paste(setdiff(random_variables, "condition"), collapse = ", "))
  }
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
#' @importFrom limma eBayes topTable lmFit
#' @export
de_ebayes = function(x, random_variables) {
  start_time <- Sys.time()

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
  if(length(rownames(random_variables)) == ncol(x) && !all(rownames(random_variables) == as.character(1:nrow(random_variables))) && rownames(random_variables) != colnames(x)) {
    append_log("input matrix x and random_variables seem misaligned; column names of x versus rownames of variables do not match", type = "error")
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
  # mask_sample_groups = match(mask_sample_groups, unique(mask_sample_groups)) - 1
  # contr_design = stats::model.matrix(~mask_sample_groups)

  fit = suppressWarnings(limma::eBayes(limma::lmFit(x, contr_design)))
  # !! sort.by="none" keeps the output table aligned with input matrix
  # always extract the first regression variable (intercept has index 1, so the column number we need is 2)
  result = suppressMessages(limma::topTable(fit, number = nrow(x), coef = 2, adjust.method = "fdr", sort.by = "none", confint = TRUE))
  ## fit$coefficients and fit$stdev.unscaled contain the intercept, remove from ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
  # note;
  result$effectsize = fit$coefficients[,2] / sqrt(fit$s2.post)
  # eBayes standard error, according to Gordon Smyth   https://support.bioconductor.org/p/70175/
  result$standarderror = sqrt(fit$s2.post) * fit$stdev.unscaled[,2]
  result$standarddeviation = sqrt(fit$s2.post)

  append_log_timestamp("eBayes", start_time)
  return(as_tibble(result) %>%
           mutate(protein_id = rownames(x)) %>%
           select(protein_id, pvalue = P.Value, qvalue = adj.P.Val, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation, standarderror))
}


#' Wrapper function for msEmpiRe::de.ana()
#'
#' input data must NOT be log transformed
#'
#' TODO: proper standarderror computation
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

  # bugfix: fix the random seed to ensure MS-EmpiRe release 1 results are the same (given the exact same input data)
  # can verify the problem case (for release 1) by disabling this line and running the respective unit-test
  set.seed(123)

  start_time <- Sys.time()
  # ms-empire expects non-logtransformed data without NA values
  x = Biobase::exprs(eset)
  if (input_intensities_are_log2) {
    x = 2^x
  }
  x[!is.finite(x)] = 0
  Biobase::exprs(eset) = x

  # protein identifiers are expected in this column
  Biobase::fData(eset)["prot.id"] = Biobase::fData(eset)[, "protein_id"]
  capture.output(result <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(eset))))
  append_log_timestamp("MS-EmpiRe", start_time)

  return(tibble(
    protein_id = as.character(result$prot.id),
    # pvalue = as.numeric(result$prot.p.val),
    # qvalue = as.numeric(result$prot.p.adj),
    pvalue = as.numeric(result$p.val),
    qvalue = as.numeric(result$p.adj),
    foldchange.log2 = as.numeric(result$log2FC),
    effectsize = as.numeric(result$log2FC) / as.numeric(result$prot.sd),
    tstatistic = NA,
    standarddeviation = as.numeric(result$prot.sd),
    standarderror = NA
  ))
}



#' Wrapper function for msqrob, using the implementation from the msqrobsum package
#'
#' @param eset a Biobase ExpressionSet that contains the peptide intensities. required attributes; fData(): protein_id and pData(): condition
#' @param use_peptide_model if true, apply msqrob. if false, apply msqrobsum
#' @param input_intensities_are_log2 whether the provided intensities are already log2 transformed
#' @param protein_rollup_robust for msqrobsum analyses, whether to use 'robust' protein rollup (TRUE) or traditional peptide intensity summation (FALSE)
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms
#'
#' @importFrom Biobase exprs pData fData
#' @importFrom MSnbase as.MSnSet.ExpressionSet
#' @export
de_msqrobsum_msqrob = function(eset, use_peptide_model = T, input_intensities_are_log2 = T, protein_rollup_robust = T, random_variables = NULL) {
  start_time <- Sys.time()

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

  msnset = MSnbase::as.MSnSet.ExpressionSet(eset)
  if (use_peptide_model) {
    form_string = c("expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id)", "expression ~ (1 | condition)")
    # form_string = "expression ~ (1 | peptide_id) + (1 | sample_id) + (1 | condition)"
  } else {
    form_string = "expression ~ (1 | condition)"
  }

  # compose updated formula
  if(length(random_variables) > 0) {
    form_string = sapply(form_string, function(x) paste(c(x, sprintf("(1 | %s)", random_variables)), collapse=" + "), USE.NAMES = FALSE)
  }

  form_string = unique(c(form_string, "expression ~ (1 | condition)"))
  append_log(paste("MSqRob regression formulas (prioritized, if model fit fails due to lack of data the next formula is tried);", paste(form_string, collapse = "  ,  ")), type = "info")

  # always add y~(1|condition) without any optional random_variables as a last-resort model if all other model specifications fail (due to lack of data) -->> then convert from string to formula/expression
  form_eval = sapply(form_string, function(x) eval(parse(text=x)), USE.NAMES = FALSE)

  # MSqRob model (re-implemented in msqrobsum package, by original authors)
  if (use_peptide_model) {
    result = suppressWarnings(msqrobsum(data = msnset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrob"))
    # result = suppressWarnings(msqrobsum(data = msnset, formulas = c(eval(parse(text=form_string)), expression ~ (1 | condition)), group_vars = "protein_id", contrasts = "condition", mode = "msqrob"))
    # our version 1, without user-specified random variables; form = c(expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id), expression ~ (1 | condition))
  } else {
    # rollup to protein level, 'robust' approach appears to be an improvement over the traditional 'sum'
    protset = suppressWarnings(suppressMessages(combineFeatures(msnset, fun = ifelse(protein_rollup_robust, "robust", "sum"), groupBy = Biobase::fData(msnset)$protein_id)))
    result = suppressWarnings(msqrobsum(data = protset, formulas = form_eval, group_vars = "protein_id", contrasts = "condition", mode = "msqrobsum"))
    ## our version 1, without user-specified random variables;
    #result = suppressWarnings(msqrobsum(data = protset, expression ~ (1 | condition), contrasts = "condition", mode = "msqrobsum", group_vars = "protein_id"))
  }

  result_unpacked = as_tibble(result) %>%
    dplyr::select(protein_id, contrasts) %>%
    tidyr::unnest(cols = contrasts)


  # sum(is.finite(msqrobsum_result_contrast$qvalue) & msqrobsum_result_contrast$qvalue <= qval_signif)
  append_log_timestamp(ifelse(use_peptide_model, "MSqRob", "MSqRobSum"), start_time)

  result_unpacked$sigma_post = result_unpacked$se / result_unpacked$sigma_contrast
  result_unpacked$effectsize = result_unpacked$logFC / result_unpacked$sigma_post

  # prepare output table in our standard format
  result = as_tibble(result_unpacked) %>%
    mutate(protein_id = as.character(protein_id)) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation = sigma_post, standarderror = se)

  if(length(random_variables) > 0) {
    result = result %>% add_column(contrast_ranvars = paste(random_variables, collapse = ", "))
  }
  return(result)
}
