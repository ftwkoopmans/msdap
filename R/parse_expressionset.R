

#' Import a label-free proteomics dataset from a peptide-level Biobase ExpressionSet object
#'
#' provided mostly for compatability with results from prior bioinformatic analyses, as it is preferred to import raw data files.
#'
#' @param eset a Biobase ExpressionSet containing peptide data
#' @param column_fdata_protein_id the column name in Biobase::fData() that holds the protein_id
#' @param column_pdata_sample_group the column name in Biobase::pData() that holds the sample group/condition
#' @param acquisition_mode the type of label-free data. options: 'dda' or 'dia'
#' @param is_log2 logical indicating whether the provided intensities are already log2 transformed
#'
#' @importFrom Biobase pData fData exprs featureNames classVersion
#' @export
import_expressionset = function(eset, column_fdata_protein_id = "prot.id", column_pdata_sample_group = "condition", acquisition_mode = 'dda', is_log2 = FALSE) {
  reset_log()
  if(!acquisition_mode %in% c("dda", "dia")) {
    append_log("`acquisition_mode` parameter can only be 'dda' or 'dia'", type = "error")
  }

  # verify input is an expressionset
  if(!"ExpressionSet" %in% names(Biobase::classVersion(eset))) {
    append_log("this function requires a dataset of type: ExpressionSet", type = "error")
  }
  # verify that user requested metadata column is in the eset
  if(!column_pdata_sample_group %in% colnames(Biobase::pData(eset))) {
    cat("column names in ExpressionSet pData:", paste(colnames(Biobase::pData(eset)), collapse = ", "), "\n")
    append_log(paste("ExpressionSet pData does not contain the sample information column you requested:", column_pdata_sample_group), type = "error")
  }
  if(!column_fdata_protein_id %in% colnames(Biobase::fData(eset))) {
    cat("column names in ExpressionSet fData:", paste(colnames(Biobase::fData(eset)), collapse = ", "), "\n")
    append_log(paste("ExpressionSet fData does not contain the sample information column you requested:", column_fdata_protein_id), type = "error")
  }


  # sample metadata
  samples = Biobase::pData(eset)
  samples$sample_id = rownames(samples) # even if there is a column named sample_id in the metadata, we over-write
  samples = as_tibble(samples) %>% mutate_all(as.character) %>% rename(group = !!column_pdata_sample_group) %>% add_column(exclude=F)
  if(!"shortname" %in% colnames(samples)) {
    samples$shortname = samples$sample_id
  }

  # protein metadata; take one entry per unique protein_id
  proteins = Biobase::fData(eset)
  proteins = proteins[!duplicated(proteins[,column_fdata_protein_id]), ]
  proteins = as_tibble(proteins) %>% mutate_all(as.character) %>% rename(protein_id = !!column_fdata_protein_id)
  if(!"fasta_headers" %in% colnames(proteins)) {
    proteins$fasta_headers = proteins$protein_id
  }

  # convert expression matrix from wide- to long-format
  tib = as_tibble(Biobase::exprs(eset)) %>%
    add_column(peptide_id = Biobase::featureNames(eset), protein_id = Biobase::fData(eset)[,column_fdata_protein_id]) %>%
    pivot_longer(cols = c(-peptide_id, -protein_id), names_to = "sample_id", values_to = "intensity") %>%
    filter(is.finite(intensity) & intensity > 0)

  # peptide sequence/modified mapping. default = peptide ID
  tib$sequence_plain = tib$sequence_modified = tib$peptide_id
  # guess which metadata column contains peptide plain/modified sequence
  col_plainseq = grep("^sequence$|^seq$|(base|plain)[._-]*(sequence|seq)|(sequence|seq)[._-]*(base|plain)", colnames(Biobase::fData(eset)), ignore.case = T, value = T)
  if(length(col_plainseq) > 0) {
    tib$sequence_plain = Biobase::fData(eset)[,col_plainseq[1]]
    append_log(paste("found plain peptide sequences in ExpressionSet fData, column:", col_plainseq), type = "info")
  }
  col_modseq = grep("(modified|mod)[._-]*(sequence|seq)|(sequence|seq)[._-]*(modified|mod)", colnames(Biobase::fData(eset)), ignore.case = T, value = T)
  if(length(col_modseq) > 0) {
    tib$sequence_plain = Biobase::fData(eset)[,col_modseq[1]]
    append_log(paste("found modified peptide sequences in ExpressionSet fData, column:", col_modseq), type = "info")
  }


  ## conform to expected properties, standardized throughout input data parsers in this codebase
  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  if(!is_log2) {
    tib$intensity = log2(tib$intensity)
  }
  tib$intensity[!is.na(tib$intensity) & tib$intensity < 0] = 0 # note; we already removed zero intensity values when importing. here, we threshold extremely low values
  tib$rt = NA # no retention time info available
  tib$detect = T # cannot differentiate between MBR and MS/MS identified
  tib$isdecoy = F # no information on decoys in expressionset

  return(list(samples=samples, proteins=proteins, peptides=tib, acquisition_mode=acquisition_mode))
}
