
#' Import a label-free proteomics dataset from OpenSWATH
#'
#' @param filename the full file path of the input file
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_openswath = function(filename, confidence_threshold = 0.01, return_decoys = FALSE, do_plot = TRUE) {
  append_log("reading OpenSWATH report...", type = "info")

  ds = import_dataset_in_long_format(filename,
                                 attributes_required = list(sample_id = "filename",
                                                            protein_id = "ProteinName",
                                                            sequence_modified = "FullPeptideName",
                                                            charge = "Charge",
                                                            rt = "norm_RT",
                                                            isdecoy = "decoy",
                                                            intensity = "Intensity",
                                                            confidence = "m_score"),
                                 attributes_optional = list(sequence_plain = "Sequence"),
                                 confidence_threshold = confidence_threshold,
                                 return_decoys = return_decoys,
                                 do_plot = do_plot)

  ds$acquisition_mode = "dia"
  return(ds)
}



#' Import a label-free proteomics dataset from Skyline
#'
#' @param filename the full file path of the input file
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by FOR DDA ONLY (acquisition_mode='dda'): if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_skyline = function(filename, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified", return_decoys = FALSE, do_plot = TRUE) {
  append_log("reading Skyline report...", type = "info")

  ds = import_dataset_in_long_format(filename,
                                       attributes_required = list(sample_id = "FileName",
                                                                  protein_id = "ProteinName",
                                                                  sequence_modified = "ModifiedSequence",
                                                                  charge = "PrecursorCharge",
                                                                  rt = "BestRetentionTime",
                                                                  intensity = "TotalArea",
                                                                  confidence = "annotation_QValue"),
                                       attributes_optional = list(isdecoy = "IsDecoy"),
                                       select_unique_precursor_per_modseq = (acquisition_mode=="dia"),
                                       confidence_threshold = confidence_threshold,
                                     return_decoys = return_decoys,
                                     do_plot = do_plot)

  ds$acquisition_mode = acquisition_mode

  if(acquisition_mode=="dda") {
    # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
    if(collapse_peptide_by == "") {
      # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
      ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = "peptide_id")
      append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
    } else {
      ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
    }
  }

  return(ds)
}



#' Import a label-free proteomics dataset from DIA-NN
#'
#' @param filename the full file path of the input file
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param use_normalized_intensities use the abundance values as-is (recommended) or those normalized by DIA-NN
#' @export
import_dataset_diann = function(filename, confidence_threshold = 0.01, use_normalized_intensities = FALSE, do_plot = TRUE) {
  append_log("reading DIA-NN report...", type = "info")
  ds = import_dataset_in_long_format(filename,
                                       attributes_required = list(sample_id = "File.Name",
                                                                  protein_id = "Protein.Group",
                                                                  sequence_modified = "Modified.Sequence",
                                                                  sequence_plain = "Stripped.Sequence",
                                                                  charge = "Precursor.Charge",
                                                                  rt = "iRT",
                                                                  # isdecoy = "",
                                                                  intensity = ifelse(use_normalized_intensities, "Precursor.Normalised", "Precursor.Quantity"),
                                                                  confidence = "Q.Value"),
                                       # attributes_optional = list(cscore = "CScore"),
                                       confidence_threshold = confidence_threshold,
                                     return_decoys = F,
                                     do_plot = do_plot)

  ds$acquisition_mode = "dia"
  return(ds)
}



#' Import a label-free proteomics dataset from Spectronaut
#'
#' No particular settings are needed, a report with peptide abundances should be made that contains these columns: R.FileName, PG.ProteinGroups, EG.IsDecoy, EG.Library, EG.Qvalue, EG.iRTEmpirical, EG.Cscore, FG.Id, FG.MS2Quantity, FG.MS2RawQuantity
#'
#'
#'
#' @param filename the full file path of the input file
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param use_normalized_intensities use the abundance values as-is (recommended) or those normalized by Spectronaut
#' @param use_irt logical indicating whether to use standardized (IRT, EG.iRTEmpirical) or the empirical (EG.RTEmpirical) retention times
#' @param remove_shared_spectronaut_proteingroups if the peptide-to-proteingroup assignments from Spectronaut are used as is (eg; you're not mapping the peptides to some spectral library downstream), you can remove peptides that map to multiple proteingroups
#'
#' @importFrom data.table chmatch
#' @export
import_dataset_spectronaut = function(filename, confidence_threshold = 0.01, use_normalized_intensities = FALSE, use_irt = TRUE,
                                      return_decoys = FALSE, remove_shared_spectronaut_proteingroups = FALSE, do_plot = TRUE) {
  append_log("reading Spectronaut report...", type = "info")

  ds = import_dataset_in_long_format(filename,
                                       attributes_required = list(sample_id = "R.FileName",
                                                                  protein_id = c("PG.ProteinAccessions", "PG.ProteinId", "EG.ProteinId", "FG.ProteinId", "PG.ProteinGroups"),
                                                                  sequence_modified = c("EG.ModifiedSequence", "FG.Id"),
                                                                  charge = c("FG.Charge", "FG.Id"),
                                                                  rt = ifelse(use_irt, "EG.iRTEmpirical", "EG.RTEmpirical"),
                                                                  isdecoy = "EG.IsDecoy",
                                                                  intensity = c("FG.MS2PeakArea", "FG.TotalPeakArea", "FG.MS2RawQuantity", "FG.Quantity"),
                                                                  cscore = "EG.Cscore",
                                                                  confidence = "EG.Qvalue"),
                                       attributes_optional = list(spectral_library = "EG.Library",
                                                                  intensity_norm = c("FG.NormalizedMS2PeakArea", "FG.NormalizedTotalPeakArea", "FG.MS2Quantity")),
                                       confidence_threshold = confidence_threshold,
                                     return_decoys = return_decoys,
                                     do_plot = do_plot)

  if(use_normalized_intensities) {
    # if user request normalized intensties from Spectronaut report, make sure this column was successfully parsed and contains a numeric value for at least 75% of target peptides (arbitrary, but should be a robust check if data is present or not)
    if(!"intensity_norm" %in% colnames(ds$peptides) || sum(is.finite(ds$peptides$intensity_norm[!ds$peptides$isdecoy])) < 0.75 * sum(!ds$peptides$isdecoy)) {
      append_log("requested normalized peptide abundances not available in spectronaut report", type = "warning")
    } else {
      ds$peptides$intensity = ds$peptides$intensity_norm
    }
  }
  # remove normalized intensity column, whether we used it or not, to prevent column name conflicts downstream
  ds$peptides$intensity_norm = NULL

  ## obsolete, as we now check for unique sequence_modified in upstream code
  # # if multiple spectral libraries were used, select 'best' match for each precursor using the overall cscore
  # is_multi_library = "spectral_library" %in% colnames(tib) && length(unique(tib$spectral_library)) > 1
  # if(is_multi_library) {
  #   append_log("multiple spectral libraries are present in the Spectronaut report, using best matching library for each precursor (measured by overall cscore).", type = "info")
  #   x = tib %>%
  #     group_by(peptide_id) %>%
  #     summarise(cscore=sum(cscore)) %>%
  #     arrange(desc(cscore))
  #   # strip spectral library from peptide_id to obtain modifiedsequence_charge (ID for precursor)
  #   x$precursor_id = gsub("##.*", "", x$peptide_id)
  #   # since we ordered from 'best' to worst', we can want to use the first occurrence of each precursor
  #   x = x[ ! duplicated(x$precursor_id), ]
  #   tib = tib %>% filter(peptide_id %in% x$peptide_id)
  # }

  ### repair uniprot protein IDs, common issue in Spectronaut reports
  # if we find a "sp " or "tr " in the protein_id column, Spectronaut replaced pipe symbols in uniprot IDs with spaces and used a pipe symbol for separating multiple proteingroups @ shared peptides
  if(any(grepl("(sp|tr) ", ds$proteins$protein_id, ignore.case=T))) {
    append_log("Spectronaut report contains malformed protein IDs (pipe symbol replaced by whitespace), fixing now...", type = "info")
    ## 1) repair protein IDs
    # if we don't want to remove shared peptides, simply remove protein-group distinction thereby effectively creating separate proteinGroups for shared peptides
    pid = ds$proteins$protein_id
    pid = gsub("|", ifelse(remove_shared_spectronaut_proteingroups, "#", ";"), pid, fixed=T)
    pid = gsub(" ", "|", pid, fixed=T)
    ## 2) update protein IDs in both protein and peptide tables
    i = data.table::chmatch(ds$peptides$protein_id, ds$proteins$protein_id)
    ds$peptides$protein_id = pid[i]
    ds$proteins$protein_id = pid
    if(remove_shared_spectronaut_proteingroups) {
      # not canonical dplyr style, but fast because we can re-use the indexing @ i
      pid_remove = grepl("#", pid, fixed=T)
      ds$proteins = ds$proteins[!pid_remove, ]
      ds$peptides = ds$peptides[!pid_remove[i], ]
      append_log("Removed all shared peptides from Spectronaut report", type = "info")
    }
  }

  ds$acquisition_mode = "dia"
  return(ds)
}



#' Import a label-free proteomics dataset from FragPipe
#'
#' @param filename the full file path of the input file
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_fragpipe = function(filename, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified", return_decoys = FALSE, do_plot = TRUE) {
  myfragpiperepair = function(DT) {
    # ### debug
    # print(head(DT))
    # hist(DT$confidence, breaks=50)
    # hist(log10(DT$confidence), breaks=50)
    # print(sum(DT$confidence >= 0.99))
    # print(sum(DT$confidence >= 0.99) / nrow(DT))
    # print(sum(DT$confidence >= 0.95))
    # print(sum(DT$confidence >= 0.95) / nrow(DT))

    ### specifically for fragpipe; strip full path + "interact-" + ".pep.xml" extension
    # example entry in input table; C:\DATA\fragpipe_test\interact-a05191.pep.xml
    DT[ , sample_id := gsub("(.*(\\\\|/)interact\\-)|(^interact\\-)|(\\.pep\\.xml$)", "", sample_id, ignore.case=F), by=sample_id]

    ### FragPipe output has "PeptideProphet Probability", which is 0~1 ranged score with 1 = 100% confident.
    # In this R package we use confidence score = qvalue (or pvalue). So as a band-aid solution we use 1-probability to approximate
    DT$confidence = 1 - DT$confidence

    ### FragPipe output lacks peptide sequences in modified sequence column if a peptiude has no modifications
    rows = DT$sequence_modified == ""
    DT$sequence_modified[rows] = DT$sequence_plain[rows]
    DT$detect = DT$confidence <= confidence_threshold

    ### summarize peptide_id*sample*id by taking the PSM with highest abundance
    data.table::setkey(DT, sequence_modified, charge, sample_id)
    # sort such that the most abundant entry for each unique key comes out on top
    data.table::setorder(DT, -intensity, na.last = TRUE)
    # simply select first entry
    DT = DT[ , index := seq_len(.N), by = .(sequence_modified, charge, sample_id)][index == 1]
    DT$index = NULL

    return(DT)
  }


  append_log("reading FragPipe PSM report...", type = "info")
  ds = import_dataset_in_long_format(filename,
                                     attributes_required = list(sample_id = "Spectrum.File",
                                                                protein_id = "Protein.ID",
                                                                sequence_modified = "Modified Peptide",
                                                                sequence_plain = "Peptide",
                                                                charge = "Charge",
                                                                rt = "Retention",
                                                                intensity = "Intensity",
                                                                confidence = "PeptideProphet.Probability"),
                                     select_unique_precursor_per_modseq = FALSE,
                                     remove_peptides_never_above_confidence_threshold = TRUE,
                                     confidence_threshold = confidence_threshold,
                                     custom_func_manipulate_DT_onload = myfragpiperepair,
                                     return_decoys = return_decoys,
                                     do_plot = do_plot)

  ds$acquisition_mode = "dda"

  # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  return(ds)
}



#' Update the peptide-to-protein mappings using a MaxQuant search result
#'
#' Simple wrapper function for user convenience; read MaxQuant proteinsgroups using import_maxquant_proteingroups(), then update the peptide-to-protein mappings in the peptide data table update_protein_mapping()
#'
#' @param dataset your dataset
#' @param path_maxquant full path to MaxQuant output ('txt') folder
#' @param remove_shared remove all shared peptides that could not be assigned to a single proteingroup as a unique- or razor-peptide ?
#' @export
update_protein_mapping_from_maxquant = function(dataset, path_maxquant, remove_shared = TRUE) {
  # TODO: input validation
  dataset$peptides = update_protein_mapping(dataset$peptides, pep2prot = import_maxquant_proteingroups(path = path_maxquant, remove_shared = remove_shared))
  dataset$proteins = dataset$proteins %>% filter(protein_id %in% dataset$peptides$protein_id)
  return(dataset)
}



#' For advanced users: update a peptide tibble with new peptide-to-protein mappings
#'
#' @param tib peptides tibble, e.g.; dataset$peptides
#' @param pep2prot a peptide-to-protein mapping table, for instance one obtained from maxquant using import_maxquant_proteingroups()
#'
#' @importFrom data.table chmatch
#' @export
update_protein_mapping = function(tib, pep2prot) {
  ### sanity check on input parameters
  if (!is.data.frame(tib) || !all(c("protein_id", "sequence_plain") %in% colnames(tib))) {
    append_log("input peptides must be a data.frame (or tibble) and contain columns protein_id and sequence_plain", type = "error")
  }
  if (!is.data.frame(pep2prot) || !all(c("protein_id", "sequence_plain") %in% colnames(pep2prot))) {
    append_log("input pep2prot mapping must be a data.frame (or tibble) and contain columns protein_id and sequence_plain", type = "error")
  }

  count_target_plainseq = length(unique(tib$sequence_plain[!tib$isdecoy]))
  # rbenchmark::benchmark(base={count_target_plainseq = length(unique(peptides$sequence_plain[!peptides$isdecoy]))},
  #                       datatable={count_target_plainseq = data.table::uniqueN(peptides$sequence_plain[!peptides$isdecoy])},
  #                       dplyr1={x=peptides %>% filter(!isdecoy) %>% summarize(cnt = n_distinct(sequence_plain))},
  #                       dplyr2={x=n_distinct(peptides %>% filter(!isdecoy) %>% pull(sequence_plain))},
  #                       dplyr3={x=peptides %>% filter(!isdecoy) %>% distinct(sequence_plain) %>% nrow()})

  ### actual mapping code
  # match only the target peptides to new protein_id mapping table 'pep2prot'
  i = data.table::chmatch(tib$sequence_plain[!tib$isdecoy], pep2prot$sequence_plain)
  # rbenchmark::benchmark(base={  i = match(tib$sequence_plain[!tib$isdecoy], pep2prot$sequence_plain)}, dt={  i = data.table::chmatch(tib$sequence_plain[!tib$isdecoy], pep2prot$sequence_plain)})
  tib$protein_id[!tib$isdecoy] = pep2prot$protein_id[i]
  # remove unmapped peptides
  tib = tib %>% filter(isdecoy | (!is.na(protein_id) & protein_id != ""))

  ### log results
  count_target_plainseq_post_filter = length(unique(tib$sequence_plain[!tib$isdecoy]))
  if(count_target_plainseq_post_filter < count_target_plainseq) {
    append_log(sprintf("%d / %d unique target (plain)sequences were not in the provided spectral library and thus removed", count_target_plainseq - count_target_plainseq_post_filter, count_target_plainseq), type = "info")
  }

  return(tib)
}



#' For advanced users: A generic function for importing data in long-format, typically used for DIA data.
#'
#' Refer to the implementation of \code{import_dataset_spectronaut} for an example on how to use this function.
#' don't forget to set either "dia" or "dda" for dataset$acquisition_mode downstream if you implement this.
#'
#' @param filename the full file path of the input file
#' @param attributes_required a list of attributes that must be present in the input file (eg; sample names, peptide IDs, intensity values, etc.)
#' @param attributes_optional a list of optional attributes that we can do without downstream, like m/z
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param remove_peptides_never_above_confidence_threshold set to true to remove all target peptides that never pass the confidence score threshold
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param select_unique_precursor_per_modseq if multiple precursors are available for a modified sequence (eg; charge 2 and 3), should we make an effort to select the best or just return all of them ?
#' @param custom_func_manipulate_DT_onload !expert use only! (there is intricate interplay between such pre-filtering and the code in this function); provide a function that manipulates the data.table after load. example implementation in the fragpipe parser.
#'
#' @importFrom data.table fread chmatch setkey setorder
#' @export
import_dataset_in_long_format = function(filename, attributes_required, attributes_optional = list(), confidence_threshold = 0.01,
                                         remove_peptides_never_above_confidence_threshold = TRUE, select_unique_precursor_per_modseq = TRUE,
                                         custom_func_manipulate_DT_onload = NULL,
                                         return_decoys = TRUE, do_plot = TRUE) {
  # sanity check on input parameters
  if (length(filename) != 1 || !is(filename, "character") || nchar(filename) < 3 || !file.exists(filename)) {
    append_log(paste("'filename' parameter is not an existing file:", filename), type = "error")
  }
  if(length(attributes_required) == 0 || !is(attributes_required, "list")) {
    append_log("'attributes_required' parameter must be a list", type = "error")
  }
  if(!is(attributes_optional, "list")) {
    append_log("'attributes_optional' parameter must be a list", type = "error")
  }
  if(any(names(attributes_required) %in% names(attributes_optional))) {
    append_log("parameters 'attributes_required' and 'attributes_optional' must not have any overlap", type = "error")
  }

  # check for minimal set of required attributes
  minimal_set_of_required = c("sample_id", "protein_id", "sequence_modified", "charge", "rt", "intensity")
  if(!all(minimal_set_of_required %in% names(attributes_required))) {
    append_log(paste("From the minimal set of peptide attributes, the following are missing;", paste(setdiff(minimal_set_of_required, names(attributes_required)), collapse = ", ")), type = "error")
  }

  # read just the first line and find indices of expected columns (named array where value = column index at input file, name = target column name)
  headers = unlist(strsplit(readLines(filename, n = 1), "\t"))
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)

  # only read columns of interest to speed up file parsing
  col_indices = c(map_required, map_optional)
  DT = data.table::fread(filename, select = unique(col_indices), check.names = T, stringsAsFactors = F)
  colnames(DT) = names(col_indices)[!duplicated(col_indices)]
  # suppose the same column from input table is assigned to multiple attributes (eg; Spectronaut FG.Id as input for both charge and modseq)
  # j = index in 'col_indices' where we should borrow content from other column
  for(j in which(duplicated(as.numeric(col_indices)))) {
    j_source_column = names(col_indices)[col_indices == col_indices[j]][1]
    DT[,names(col_indices)[j]] = DT[[j_source_column]]
  }


  ### inject function for custom manipulation of DT. example; for some input format, repair a column. for FragPipe, the modified peptide column is empty for peptides without a modification, so there is no column with proper modseq available. have to 'manually repair', but other than that we want to maintain this generic function
  if(length(custom_func_manipulate_DT_onload) == 1 && !is.null(custom_func_manipulate_DT_onload) && is.function(custom_func_manipulate_DT_onload)) {
    DT = custom_func_manipulate_DT_onload(DT)
  }


  ### first, see what data we can remove, and perform more expensive operations after
  if("isdecoy" %in% colnames(DT)) {
    DT$isdecoy = DT$isdecoy %in% c(TRUE, 1, "decoy", "true")
  } else {
    DT$isdecoy = FALSE
  }

  # no confidence scores in the input table, assume all are 'detected'. otherwise, test confidence score vs threshold
  if(!"confidence" %in% colnames(DT)) {
    append_log("input dataset lacks peptide confidence scores, thus assuming all peptides are 'detected'", type = "warning")
    DT$detect = TRUE
  } else {
    # define detect using a threshold on confidence scores
    DT[ , detect := is.finite(confidence) & confidence <= confidence_threshold ]
  }


  # remove entries that lack crucial data
  # remove target peptides that don't have intensity value. for decoys, intensity values are not required
  DT = DT[sample_id != "" & protein_id != "" & sequence_modified  != "" & (isdecoy | (is.finite(intensity) & intensity > 0))]


  # remove path and extension from filename. using grouping, we only have to apply the regex on unique elements (making this very fast)
  DT[ , sample_id := gsub("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia)(\\.gz|\\.dia){0,1}$)", "", sample_id, ignore.case=T), by=sample_id]
  ## benchmarks show why we use data.table for all downstream string manipulations; good combination of readable code and computational performance
  # rbenchmark::benchmark(naive={s = gsub("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia)(\\.gz|\\.dia){0,1}$)", "", tib$sample_id, ignore.case=T)},
  #                       hardcoded={
  #                         sid = unique(DT$sample_id)
  #                         sid_new = gsub("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia)(\\.gz|\\.dia){0,1}$)", "", sid, ignore.case=T)
  #                         i = data.table::chmatch(DT$sample_id, sid)
  #                       },
  #                       datatable={DT[ , sample_id := gsub("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia)(\\.gz|\\.dia){0,1}$)", "", sample_id, ignore.case=T), by=sample_id]} )


  # assume the charge is either provided as a plain number, or as trailing character in a string (eg; FragmentGroup ID @ spectronaut)
  if(typeof(DT$charge) == "character") { # if input data is not numeric... (data.table::fread() should have automatically recognise it)
    DT[ , charge := sub(".*\\D(\\d+)$", "\\1", charge), by=charge]
  }
  DT[ , charge := replace_na(suppressWarnings(as.integer(charge)), 0), by=charge]

  # remove leading and trailing annotations from modified sequence, such as underscores or charge states. regex; start with non-character or opening bracking (modification), and analogous for trailing
  DT[ , sequence_modified := gsub("^_|_[^A-Z]*$", "", sequence_modified), by=sequence_modified]
  # DT[ , sequence_modified := gsub("(^[^a-zA-Z([]+)|([^a-zA-Z)\\]]+$)", "", sequence_modified, perl=T), by=sequence_modified] # universal, but a very slow regex
  # test input for regex; gsub("(^[^a-zA-Z([]+)|([^a-zA-Z)\\]]+$)", "", c("aa", "aa_", "_aa", "_(b)aa(c)_", "(b)aa(c)..", "_[b]aa[c]_", "2.[b]aa[c]", "_[c(unimod:1)]aaa"), perl=T)

  # if no plain sequence is available in input table, extract it from the modified sequence by removing everything between brackets. either [] or ()
  # note that above we validated sequence_modified is never empty. so if there the provided sequence_plain has any empty values, default back to a manipulation of modified sequences
  if (!"sequence_plain" %in% colnames(DT) || any(DT$sequence_plain == "")) {
    DT[ , sequence_plain := gsub("(\\[[^]]*\\])|(\\([^)]*\\))", "", sequence_modified), by=sequence_modified]
  }

  # peptide_id = modified sequence + charge (eg. unique precursors) + spectral library (so we can support scenario's where the raw data was matched to multiple spectral libraries)
  if("spectral_library" %in% colnames(DT)) {
    DT[ , peptide_id := paste0(sequence_modified, "_", charge, "###", spectral_library), by= .(sequence_modified, charge)]
  } else {
    DT[ , peptide_id := paste0(sequence_modified, "_", charge), by= .(sequence_modified, charge)]
  }

  # peptide_id is a key, it must be unique per sample
  if(anyDuplicated(DT[isdecoy==FALSE], by=c("peptide_id", "sample_id")) != 0) {
    append_log("error: peptide_id*sample_id combinations are not unique (eg; same peptide_id occurs multiple times per sample_id). If multiple data points per precursor*sample are available in your dataset, make sure to provide proper columns with the modified sequence, charge and spectral library that each data points originates from!", type = "error")
  }

  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  DT$intensity = log2(DT$intensity)
  DT$intensity[!is.finite(DT$intensity)] = NA
  if("intensity_norm" %in% colnames(DT)) {
    DT$intensity_norm = log2(DT$intensity_norm)
    DT$intensity_norm[!is.finite(DT$intensity_norm)] = NA
  }

  # for DIA data we can compute Cscore histograms of both target and decoy peptides (which will be shown in report downstream)
  plotlist = list()
  if(do_plot && "cscore" %in% colnames(DT)) {
    plotlist$ggplot_cscore_histograms = plot_dia_cscore_histograms(as_tibble(DT))
  }

  if(!return_decoys) {
    DT = DT[isdecoy==FALSE]
  }

  # remove all target peptides that do not have a single detect (eg; DIA dataset, some peptide from the spectral library that is never confidently identified/quantified)
  if(remove_peptides_never_above_confidence_threshold) {
    DT_valid = DT[detect==TRUE, .(peptide_id), by=peptide_id]
    DT = DT[isdecoy==TRUE | data.table::chmatch(peptide_id, DT_valid$peptide_id) != 0]
  }


  # remove peptides attributed to multiple protein identifiers
  DT_valid = DT[isdecoy==FALSE][ , .(sequence_plain), by = .(sequence_plain, protein_id)][ , .( fail = .N != 1), by = sequence_plain]
  if(any(DT_valid$fail)) {
    append_log(sprintf("%d unique target (plain)sequences ambiguously mapped to multiple proteins and thus removed. Examples; %s",
                       sum(DT_valid$fail), paste(head(DT_valid$sequence_plain[DT_valid$fail], 5), collapse=", ")), type = "info")
    # keep entries that are either decoys or their plain sequence matches the subset that did NOT fail the above uniqueness test
    DT = DT[isdecoy==TRUE | data.table::chmatch(sequence_plain, (DT_valid[fail==FALSE])$sequence_plain) != 0]
  }

  # select 'best' peptide_id for each sequence_modified, in case there are multiple data points provided due to multiple charges, spectral libraries, etc.
  if(select_unique_precursor_per_modseq) {
    # select 'best' peptide_id for each sequence_modified, in case there are multiple data points provided due to multiple charges, spectral libraries, etc.
    data.table::setkey(DT, peptide_id)
    # simple summary statistics per peptide_id
    x = DT[isdecoy==FALSE][ , .(sequence_modified=sequence_modified[1], ndetect=sum(detect==TRUE), abundance=mean(intensity)), by=peptide_id]
    nprecursor = nrow(x)
    # sort such that the desireable sequence_modified for each peptide_id comes out on top
    data.table::setorder(x, -ndetect, -abundance, na.last = TRUE)
    # simply select first entry
    x = x[ , index := seq_len(.N), by = sequence_modified][index == 1]
    nmodseq = nrow(x)
    # debug; x[sequence_modified=="KNPDSQYGELIEK"]; y[sequence_modified=="KNPDSQYGELIEK"]

    DT = DT[isdecoy==TRUE | data.table::chmatch(peptide_id, x$peptide_id) != 0]
    if(nprecursor != nmodseq) {
      append_log(sprintf("%d/%d precursors remain after selecting the 'best' precursor for each modified sequence", nmodseq, nprecursor), type = "info")
    }
  }


  if(!any(DT$isdecoy == FALSE)) {
    append_log("there must be target peptides ('isdecoy' column == FALSE) in the peptides tibble", type = "error")
  }

  proteins = empty_protein_tibble(DT)

  # TODO: log_peptide_tibble_pep_prot_counts(tib)
  return(list(peptides=tibble_peptides_reorder(as_tibble(DT)), proteins=proteins, plots=plotlist))
}
