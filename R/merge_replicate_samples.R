
#' Merge technical replicates prior to downstream analysis
#'
#' @description
#' Replicate measurements of the same "biological sample ID" are merged by this function, updating the peptide and
#' sample tables; foreach set of rawfiles (sample_id in the sample metadata table) from the same biological sample,
#' all but 1 sample_id will be removed and the remaining sample_id will be assigned the mean peptide values
#' (intensity, retention time, etc.) across all replicates.
#'
#' Peptide log2 intensity values can be rescaled prior to merging values across samples to prevent potential bias
#' towards samples with higher sample loading (recommended, default setting).
#'
#' Importantly, the sample metadata table must contain a column with unique identifiers for all samples that must be merged (so NOT replicate numbers!). Keep a close eye on the output log (also found in report.pdf) to confirm the right samples have been merged !
#'
#'
#' Below example shows a sample metadata table with 2 technical replicates for mouse1, 2 technical replicates for mouse2, and no technical replicates for the mice in group B (bioid column is left empty).
#'
#' - sample_id      group   bioid
#' - rawfile1      grpA      mouse1
#' - rawfile2      grpA      mouse1
#' - rawfile3      grpA      mouse2
#' - rawfile4      grpA      mouse2
#' - rawfile5      grpB
#' - rawfile6      grpB
#'
#' @examples \dontrun{
#'   # first, import your dataset and sample metadata
#'   dataset = import_dataset_diann(...)
#'   dataset = import_sample_metadata(...)
#'
#'   # merge technical replicates
#'   dataset = merge_replicate_samples(
#'     dataset,
#'     colname = "bioid",
#'     minsample = 2,
#'     rescale_intensities = TRUE
#'   )
#'
#'   # proceed with typical MS-DAP usage
#'   dataset = analysis_quickstart(...)
#' }
#' @param dataset dataset with sample metadata attached
#' @param colname column name in the samples table that represents identifiers for samples that belong together / should be merged; give the same sample-unique label/ID to each row in the sample metadata that should be merged. NOT A COLUMN WITH TECHNICAL REPLICATES IDS !
#' @param minsample minimum number of 'technical replicates' samples where a peptide must be observed. e.g. if set to 2 and there are 3 technical replicates, this will retain all peptides that have an intensity value in at least 2 out of 3 replicates. If you provide a value larger than the number of biological replicates, this limit will be capped at the number of available replicates (e.g. min_sample=3 for biosamples with 2 replicates -> retain peptides available in 2/2 replicates). By default, set to 1 (retaining all peptides found in any of the technical replicates)
#' @param rescale_intensities boolean value indicating whether by-sample median normalization should be applied prior to averaging log2 peptide intensities across samples. Strongly recommended if you didn't import data that was already normalized at peptide-level. Doesn't hurt much if normalized data is used as input, so this is enabled by default
#' @export
merge_replicate_samples = function(dataset, colname, minsample = 1L, rescale_intensities = TRUE) {
  if( ! (length(rescale_intensities) == 1 && rescale_intensities %in% c(TRUE, FALSE)) ) {
    append_log("merge_replicate_samples() parameter 'rescale_intensities' must be a single boolean value", type = "error")
  }
  if( ! (length(minsample) == 1 && is.finite(minsample) && is.numeric(minsample) && minsample >= 1) ) {
    append_log("merge_replicate_samples() parameter 'minsample' must be a single integer value >= 1", type = "error")
  }
  if(!is.integer(minsample)) {
    minsample = as.integer(round(minsample))
  }

  # nothing to do
  if(length(colname) == 0 || anyNA(colname)) {
    return(dataset)
  }
  # unknown column
  if(length(colname) != 1 || !is.character(colname) || colname == "" || !tolower(colname) %in% tolower(colnames(dataset$samples))) {
    append_log(paste0(
      "to merge technical replicate samples prior to downstream analysis, the function merge_replicate_samples() requires a column name in the sample metadata table that describes samples that belong to the same biological replicate; the provided column named '",
      colname,
      "' was not found in dataset$samples (to disable merging of raw files, provide NULL or NA for the colname parameter)"
    ), type = "error")
  }
  # guard against invalid column names / don't merge these
  if(tolower(colname) %in% c("group", "exclude")) {
    append_log("to merge technical replicate samples prior to downstream analysis, provide a column name in the sample metadata table that describes samples that belong to the same biological replicate. Provided column name is illegal (cannot be 'group' or 'exclude')", type = "error")
  }

  # case insensitive matching (works because we already checked for presence of 'colname' upstream)
  if( ! colname %in% colnames(dataset$samples)) {
    colname = colnames(dataset$samples)[match(tolower(colname), tolower(colnames(dataset$samples)))]
  }

  # collect sample_id for which a biological sample ID is provided
  sample_ref = dataset$samples %>%
    select(sample_id, group, samesample = !!colname) %>%
    mutate(samesample = as.character(samesample)) %>%
    filter(!is.na(samesample) & samesample != "")
  # remove samples that have a unique biological sample ID, i.e. there are no technical replicates
  sample_ref = sample_ref %>%
    filter(samesample %in% (sample_ref %>% count(samesample) %>% filter(n > 1) %>% pull(samesample)) )

  # nothing to do
  if(nrow(sample_ref) == 0) {
    append_log(paste0(
      "no technical replicates found; all biological sample identifiers that are defined in the '", colname,
      "' column of the sample metadata table are unique so there is nothing to do.\nDid you mistakenly enter 'replicate numbers' instead of an ID per biological sample? Provide the same identifier across all (technical replicate) samples that need to be merged."
    ), type = "warning")
    return(dataset)
  }

  # samesample IDs are spread across groups, user must have made an error in input
  sample_ref_grouperror = sample_ref %>% distinct(group, samesample) %>% count(samesample) %>% filter(n > 1)
  if(nrow(sample_ref_grouperror) > 0) {
    append_log(paste0(
      "invalid definition of 'biological sample identifiers' is used as input for merging technical replicates; column '",
      colname,
      "' contains the same biological sample identifier across multiple sample groups (samples that are technical replicates must not have different values in the 'group' column)\nDid you mistakenly enter 'replicate numbers' instead of an ID per biological sample? Provide the same identifier across all (technical replicate) samples that need to be merged."
    ), type = "error")
  }



  # all duplicated 'samesample' entries should be merged
  biorep_ids = unique(sample_ref$samesample)
  for(biorep_id in biorep_ids) {
    # length of sid is always at least 2; we removed unique 'samesample' upstream
    sid = sample_ref %>% filter(samesample == biorep_id) %>% pull(sample_id)
    # rescale sample up-front to ensure that missing values aren't biased towards the sample that doesn't have the missing value
    # e.g. sample A is more abundant in input than sample B.  merging intensity values and assuming MCAR for NA values, a value missing in B will have inflated abundance as compared to a peptides without missing values (mean over both samples)
    pep_subset = dataset$peptides %>% filter(sample_id %in% sid & is.finite(intensity))
    # peptide*sample matrix
    m_int = pep_subset %>%
      select(peptide_id, sample_id, intensity) %>%
      tidyr::pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "intensity") %>%
      as_matrix_except_first_column()

    # instead of scaling samples such that median is zero, adjust by mean shift so output values are of the same order as input
    # example data for verifying the rescaling code; m_int=cbind(1:6, 2:7, 0:5)
    if(rescale_intensities) {
      sid_scale = matrixStats::colMedians(m_int, na.rm = TRUE)
      sid_scale = sid_scale - mean(sid_scale)
      for(j in 1:ncol(m_int)) {
        m_int[,j] = m_int[,j] - sid_scale[j]
      }
    }

    # new peptide table
    i = match(rownames(m_int), pep_subset$peptide_id)
    tib_new = tibble::tibble(
      sample_id = sid[1],
      peptide_id = rownames(m_int),
      protein_id = pep_subset$protein_id[i],
      sequence_plain = pep_subset$sequence_plain[i],
      sequence_modified = pep_subset$sequence_modified[i],
      intensity = matrixStats::rowMeans2(m_int, na.rm = TRUE)
    )

    # remove peptides not observed in enough replicates by setting value to NA
    # (don't subset the table, downstream code relies on preservation of table order)
    if(minsample > 1) {
      # importantly, threshold filtering to the number of available replicates
      rows = matrixStats::rowSums2(is.finite(m_int)) < min(minsample, ncol(m_int))
      tib_new$intensity[rows] = NA
      if(all(rows)) {
        append_log(sprintf(
          "when merging %d technical replicates for biological sample ID '%s' (described in column '%s' of the sample metadata table), the 'minsample=%d' filtering setting resulted in zero peptides (i.e. all data for this biosample was removed). Consider setting a lower value for 'minsample'",
          length(sid), biorep_id, colname, minsample), type = "warning")
      }
    }


    # note that the pivot_wider is always on the same input tibble, thus the order of output rownames is always the same as above code for intensities
    classes = sapply(pep_subset, typeof)
    # skip all intensity columns  +  unsupported types   (the former also ensures the earlier minsample filter is not overwritten)
    classes = classes[ ! grepl("^intensity", names(classes), ignore.case = TRUE) & classes %in% c("logical", "integer", "double", "numeric")] # don't use intersect() as it drops the array names
    mybooltest = function(x) { is.finite(x) & x > 0 }
    for(i in seq_along(classes)) {
      i_col = names(classes)[i]
      i_type = unname(classes[i])
      if(i_type == "logical") {
        tib_new[,i_col] = pep_subset %>%
          select(peptide_id, sample_id, value = !!i_col) %>%
          mutate(value = as.integer(value %in% TRUE)) %>%
          tidyr::pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "value", values_fill = list(value = 0L)) %>%
          as_matrix_except_first_column() %>%
          matrixStats::rowSums2(na.rm = TRUE) %>%
          mybooltest()
      }
      if(i_type %in% c("integer", "double", "numeric")) {
        tib_new[,i_col] = pep_subset %>%
          select(peptide_id, sample_id, value = !!i_col) %>%
          tidyr::pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "value") %>%
          as_matrix_except_first_column() %>%
          matrixStats::rowMeans2(na.rm = TRUE)
      }
    }


    tib_new = tib_new %>%
      # remove peptides without a 'merged intensity'
      filter(is.finite(intensity)) %>%
      mutate_all(unname) %>%
      # arrange columns in same order as reference
      select(order(match(names(.), names(dataset$peptides))))

    # all samples unrelated to current set of samples-to-merge  +  new data
    # note that tib_new can be empty (e.g. no peptides present in both replicates @ min_sample=2) --> we still want to execute this line to remove original sample_id
    dataset$peptides = bind_rows(dataset$peptides %>% filter( ! sample_id %in% sid), tib_new)


    append_log(sprintf("merged %d technical replicates for biological sample ID '%s' (described in column '%s' of the sample metadata table);\n%s",
                       length(sid), biorep_id, colname, paste(sid, collapse = "\n")), type = "info")
  }




  ### finally, update the sample table; remove all sample_id for technical replicates that aren't used a "representative sample_id" per biosample
  # count number of times each biorep is observed into a table with columns 'samesample' and 'n'
  tmp = sample_ref %>% count(samesample)
  # align with samples table, then replace missing values with 1 (i.e. all samples not specified for merging will be listed as '1 replicate')
  tmp2 = tmp$n[match(dataset$samples %>% select(!!colname) %>% pull(), tmp$samesample)]
  tmp2[!is.finite(tmp2) | tmp2 < 2] = 1L
  colname_repcount = paste0(colname, "_repcount")
  dataset$samples[,colname_repcount] = as.integer(tmp2)

  dataset$samples = dataset$samples %>%
    ungroup() %>%
    # enforce character type
    mutate_at(dplyr::vars(colname), as.character) %>%
    # move the column name that holds the counts next to the 'bioid' definition
    relocate(!!colname_repcount, .after = colname) %>%
    # remove obsoleted sample_id
    filter(sample_id %in% unique(dataset$peptides$sample_id))

  # reset sample indices and update the counts per sample
  dataset$samples = dataset$samples %>% mutate(sample_index = 1L:nrow(dataset$samples))
  dataset$samples = peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))

  return(dataset)
}
