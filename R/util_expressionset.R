
#' placeholder title
#' @param peptides todo
#' @param proteins todo
#' @param samples todo
#'
#' @importFrom Biobase pData fData exprs annotatedDataFrameFrom
tibble_as_eset = function(peptides, proteins, samples) {
  # append_log("grouping peptides in an ExpressionSet...", type = "info")

  # from long format to wide; peptide intensity matrix
  peptide_intensities = as_matrix_except_first_column(peptides %>%
                                                        select(peptide_id, sample_id, intensity) %>%
                                                        spread(sample_id, intensity))

  # order the columns/samples as specified in `samples`
  j = order(match(colnames(peptide_intensities), samples$sample_id))  # example; order(match(c("a", "c", "x", "b"), letters[1:3]))
  peptide_intensities = peptide_intensities[,j,drop=F]

  # create ExpressionSet with empty feature/protocol data
  eset = Biobase::ExpressionSet(
    assayData = peptide_intensities,
    featureData = Biobase::annotatedDataFrameFrom(peptide_intensities, byrow = T),
    protocolData = Biobase::annotatedDataFrameFrom(peptide_intensities, byrow = F)
  )

  ### sample metadata
  # align the sample order between the eset and the source data
  sid = rownames(Biobase::pData(eset))
  Biobase::pData(eset) = data.frame(samples %>% slice(match(sid, sample_id)), row.names = sid)

  ### peptides -->> protein metadata
  pid = rownames(Biobase::fData(eset))
  pep = peptides %>%
    select(peptide_id, sequence_plain, sequence_modified, protein_id) %>%
    slice(match(pid, peptide_id)) %>%
    left_join(proteins, by = "protein_id")

  Biobase::fData(eset) = data.frame(pep, row.names = pid)
  # Biobase::fData(eset) = data.frame(lapply(pep, as.factor), row.names = pid) # if downstream code expects factors in ExpressionSet

  return(eset)
}



#' Convert a peptide-level ExpressionSet to protein-level by summing their respective peptide intensities per sample
#'
#' @param eset_peptides A Biobase ExpressionSet with peptide data. Intensity values are assumed to be log2 transformed! fData() must contain protein_id attributes
#' @importFrom Biobase fData pData exprs
#' @export
eset_from_peptides_to_proteins = function(eset_peptides) {
  fdata = Biobase::fData(eset_peptides)
  if(!("protein_id" %in% colnames(fdata))) {
    append_log("expressionset fData() must contain column 'protein_id'", type = "error")
  }

  # peptide abundance matrix -->> row-wise aggregation. don't replace NA with zero before undoing the log2 transformation
  expr_prot = as_tibble(2^Biobase::exprs(eset_peptides)) %>% replace(is.na(.), 0) %>% add_column(protein_id=fdata$protein_id) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA)
  m = as.matrix(expr_prot[,-1])
  rownames(m) = expr_prot$protein_id
  # convert to log2 again (note; we already replaced zero's with NA above)
  m = log2(m)

  # create ExpressionSet with empty feature/protocol data
  eset_proteins = Biobase::ExpressionSet(
    assayData = m,
    featureData = Biobase::annotatedDataFrameFrom(m, byrow = T),
    protocolData = Biobase::annotatedDataFrameFrom(m, byrow = F)
  )

  Biobase::pData(eset_proteins) = Biobase::pData(eset_peptides)
  Biobase::fData(eset_proteins) = fdata[match(expr_prot$protein_id, fdata$protein_id), !grepl("sequence|peptide", colnames(fdata), ignore.case = T)]

  ## reference implementation, v1
  #' @importFrom MSnbase as.MSnSet.ExpressionSet combineFeatures
  # # need non-log values for combineFeatures()
  # tmp_eset = eset_peptides
  # Biobase::exprs(tmp_eset) = 2^Biobase::exprs(tmp_eset)
  # Biobase::exprs(tmp_eset)[!is.finite(Biobase::exprs(tmp_eset))] = 0
  #
  # # rollup to protein level, for eBayes we use the traditional sum approach
  # eset_proteins = suppressWarnings(suppressMessages(MSnbase::combineFeatures(MSnbase::as.MSnSet.ExpressionSet(tmp_eset), fun = "sum", groupBy = Biobase::fData(tmp_eset)$protein_id)))
  #
  # # convert to log2 again
  # Biobase::exprs(eset_proteins) = log2(Biobase::exprs(eset_proteins))
  # Biobase::exprs(eset_proteins)[!is.finite(Biobase::exprs(eset_proteins))] = NA

  return(eset_proteins)
}



#' placeholder title
#' @param eset todo
#' @param valid_peptide_ids todo
#' @param valid_sample_ids todo
#'
#' @importFrom Biobase fData sampleNames
subset_eset = function(eset, valid_peptide_ids, valid_sample_ids) {
  fdata = Biobase::fData(eset)
  if(!("peptide_id" %in% colnames(fdata))) {
    append_log("expressionset fData() must contain column 'peptide_id'", type = "error")
  }

  cols = rep(T, ncol(fdata))
  if (length(valid_sample_ids) > 0 && all(!is.na(valid_sample_ids) && valid_sample_ids != "")) {
    cols = Biobase::sampleNames(eset) %in% valid_sample_ids
  }
  # simply use matrix-style subsetting.  documentation @ https://genomicsclass.github.io/book/pages/eset.html
  return(eset[fdata$peptide_id %in% valid_peptide_ids, cols])
}
