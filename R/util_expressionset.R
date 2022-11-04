
#' peptide-level ExpressionSet from long-format peptide data
#'
#' @param peptides peptide tibble
#' @param proteins protein tibble
#' @param samples sample tibble
#' @importFrom Biobase pData fData exprs annotatedDataFrameFrom
#' @export
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



#' protein-level ExpressionSet from protein abundance matrix
#'
#' @param x protein abundance matrix
#' @param eset A Biobase ExpressionSet with protein_id and sample_id metadata
#' @importFrom Biobase fData pData
#' @export
protein_eset_from_data = function(x, eset) {
  fdata = Biobase::fData(eset)
  pdata = Biobase::pData(eset)
  if(!("protein_id" %in% colnames(fdata))) {
    append_log("expressionset fData() must contain column 'protein_id'", type = "error")
  }
  if(!all(rownames(x) %in% fdata$protein_id)) {
    append_log("expressionset fData() column 'protein_id' must contain an entry for each 'rowname' of matrix x", type = "error")
  }
  if(ncol(x) != nrow(pdata) || !all(colnames(x) == pdata$sample_id)) {
    append_log("expressionset pData() must align with matrix x; each row in pdata must describe a column in x AND the pdata column 'sample_id' must match column names of x", type = "error")
  }

  # create ExpressionSet with empty feature/protocol data
  eset_proteins = Biobase::ExpressionSet(
    assayData = x,
    featureData = Biobase::annotatedDataFrameFrom(x, byrow = T),
    protocolData = Biobase::annotatedDataFrameFrom(x, byrow = F)
  )

  # transfer metadata
  Biobase::fData(eset_proteins) = fdata[match(rownames(x), fdata$protein_id), !grepl("sequence|peptide", colnames(fdata), ignore.case = T), drop=F]
  Biobase::pData(eset_proteins) = pdata
  return(eset_proteins)
}

