
#' placeholder title
#' @param peptides todo
#' @param samples todo
#' @param eset_peptides todo
#' @param eset_proteins todo
#' @param input_intensities_are_log2 todo
#'
#' @importFrom Biobase exprs pData
de_interface_ebayes = function(peptides=NULL, samples=NULL, eset_peptides=NULL, eset_proteins=NULL, input_intensities_are_log2 = T) {
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

  de_ebayes(x = Biobase::exprs(eset_proteins), mask_sample_groups = Biobase::pData(eset_proteins)$condition) %>% add_column(algo_de = "ebayes")
}



#' Wrapper function for limma::ebayes()
#'
#' rownames of the data matrix (log2 intensities) must be the protein_id for downstream compatibility
#'
#' @param x log2 transformed protein intensity matrix
#' @param mask_sample_groups numeric or string array indicating to which group each column/sample belongs. NO NA VALUES !
#'
#' @importFrom stats model.matrix
#' @importFrom limma eBayes topTable
#' @export
de_ebayes = function(x, mask_sample_groups) {
  if(any(is.na(mask_sample_groups))) {
    append_log("mask_sample_groups must not contain any NA values", type = "error")
  }
  if(length(rownames(x)) == 0) {
    append_log("input matrix must have rownames (these are assumed to be protein identifiers)", type = "error")
  }

  start_time <- Sys.time()
  mask_sample_groups = match(mask_sample_groups, unique(mask_sample_groups)) - 1

  # ref implementation: contr_design = model.matrix(~ rep(0:1, c(length(contr_cols_a), length(contr_cols_b))))
  contr_design = stats::model.matrix(~mask_sample_groups)
  fit = suppressWarnings(limma::eBayes(lmFit(x, contr_design)))
  # !! sort.by="none" keeps the output table aligned with input matrix
  result = suppressMessages(limma::topTable(fit, number = nrow(x), adjust.method = "fdr", sort.by = "none", confint = TRUE))
  ## fit$coefficients and fit$stdev.unscaled contain the intercept, remove from ES and SE computation
  # eBayes effect size: Cohen's d in limma, according to Gordon Smyth   https://support.bioconductor.org/p/71747/#71781
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
#'
#' @importFrom Biobase exprs pData fData
#' @importFrom MSnbase as.MSnSet.ExpressionSet
#' @export
de_msqrobsum_msqrob = function(eset, use_peptide_model = T, input_intensities_are_log2 = T, protein_rollup_robust = T) {
  start_time <- Sys.time()

  if(!("protein_id" %in% colnames(Biobase::fData(eset)))) {
    append_log("ExpressionSet fData() must contain column 'protein_id'", type = "error")
  }
  if(use_peptide_model && !("peptide_id" %in% colnames(Biobase::fData(eset)))) {
    append_log("ExpressionSet fData() must contain column 'peptide_id' (when using peptide-level msqrob model)", type = "error")
  }
  if(!all(c("condition", "sample_id") %in% colnames(Biobase::pData(eset)))) {
    append_log("ExpressionSet pData() must contain columns 'condition' and 'sample_id'", type = "error")
  }

  # transform to log2 if input data is non-log
  if (!input_intensities_are_log2) {
    x = log2(Biobase::exprs(eset))
    x[!is.finite(x)] = NA
    Biobase::exprs(eset) = x
  }

  msnset = MSnbase::as.MSnSet.ExpressionSet(eset)

  if (use_peptide_model) {
    # MSqRob model (re-implemented in msqrobsum package, by original authors)
    form = c(expression ~ (1 | condition) + (1 | sample_id) + (1 | peptide_id), expression ~ (1 | condition))
    result = suppressWarnings(msqrobsum(data = msnset, formulas = form, group_vars = c("protein_id"), contrasts = "condition", mode = "msqrob"))
  } else {
    # rollup to protein level, 'robust' approach appears to be an improvement over the traditional 'sum'
    protset = suppressWarnings(suppressMessages(combineFeatures(msnset, fun = ifelse(protein_rollup_robust, "robust", "sum"), groupBy = Biobase::fData(msnset)$protein_id)))
    result = suppressWarnings(msqrobsum(data = protset, expression ~ (1 | condition), contrasts = "condition", mode = "msqrobsum", group_vars = "protein_id"))
  }

  result_unpacked = as_tibble(result) %>%
    dplyr::select(protein_id, contrasts) %>%
    tidyr::unnest(cols = contrasts)


  # sum(is.finite(msqrobsum_result_contrast$qvalue) & msqrobsum_result_contrast$qvalue <= qval_signif)
  append_log_timestamp(ifelse(use_peptide_model, "MSqRob", "MSqRobSum"), start_time)

  result_unpacked$sigma_post = result_unpacked$se / result_unpacked$sigma_contrast
  result_unpacked$effectsize = result_unpacked$logFC / result_unpacked$sigma_post

  # prepare output table in our standard format
  return(as_tibble(result_unpacked) %>%
    mutate(protein_id = as.character(protein_id)) %>%
    select(protein_id, pvalue, qvalue, foldchange.log2 = logFC, effectsize, tstatistic = t, standarddeviation = sigma_post, standarderror = se))
}
