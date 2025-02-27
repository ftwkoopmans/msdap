
#' Import a label-free proteomics dataset from OpenSWATH
#'
#' @param filename the full file path of the input file
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', default: 0.01 (target value must be lesser than or equals)
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_openswath = function(filename, confidence_threshold = 0.01, return_decoys = FALSE, do_plot = TRUE) {
  reset_log()
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
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', default: 0.01 (target value must be lesser than or equals)
#' @param collapse_peptide_by FOR DDA ONLY (acquisition_mode='dda'): if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_skyline = function(filename, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified", return_decoys = FALSE, do_plot = TRUE) {
  reset_log()
  append_log("reading Skyline report...", type = "info")
  stopifnot(acquisition_mode %in% c("dda","dia"))

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
#' @param filename the full file path of the DIA-NN report; this can be either the report.tsv file or the report.parquet file (new since DIA-NN 1.9.1)
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', default: 0.01 (target value must be lesser than or equals)
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param use_irt logical indicating whether to use standardized (IRT) or the empirical (RT) retention times
#' @param use_normalized_intensities use the abundance values as-is (recommended) or those normalized by DIA-NN
#' @export
import_dataset_diann = function(filename, confidence_threshold = 0.01, use_normalized_intensities = FALSE, use_irt = TRUE, do_plot = TRUE) {
  reset_log()
  append_log("reading DIA-NN report...", type = "info")

  rt_col = ifelse(use_irt, "iRT", "RT")

  # if the filename ends with ".parquet", use the Apache arrow library to read the DIA-NN report in parquet format
  f = df = NULL
  if(grepl("\\.parquet$", filename)) {
    append_log(paste("input file:", filename), type = "info")
    df = arrow::read_parquet(
      filename,
      col_select = NULL,
      as_data_frame = TRUE,
      props = arrow::ParquetArrowReaderProperties$create(),
      mmap = TRUE
    )
  } else { # not a parquet file, assume input file is a .tsv file
    f = filename
  }

  ds = import_dataset_in_long_format(filename = f,
                                     x = df,
                                     attributes_required = list(sample_id = c("File.Name", "Run"),
                                                                protein_id = "Protein.Group",
                                                                sequence_modified = "Modified.Sequence",
                                                                sequence_plain = "Stripped.Sequence",
                                                                charge = "Precursor.Charge",
                                                                rt = rt_col,
                                                                # isdecoy = "",
                                                                intensity = ifelse(use_normalized_intensities, "Precursor.Normalised", "Precursor.Quantity"),
                                                                confidence = "Q.Value",
                                                                confidence_score = "Evidence"),
                                     # attributes_optional = list(confidence_score = "Evidence"),
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
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', default: 0.01 (target value must be lesser than or equals)
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param use_normalized_intensities use the abundance values as-is (recommended) or those normalized by Spectronaut
#' @param use_irt logical indicating whether to use standardized (IRT, EG.iRTEmpirical) or the empirical (EG.RTEmpirical) retention times
#' @param remove_shared_spectronaut_proteingroups if the peptide-to-proteingroup assignments from Spectronaut are used as is (eg; you're not mapping the peptides to some spectral library downstream), you can remove peptides that map to multiple proteingroups
#'
#' @importFrom data.table chmatch "%chin%"
#' @export
import_dataset_spectronaut = function(filename, confidence_threshold = 0.01, use_normalized_intensities = FALSE, use_irt = TRUE,
                                      return_decoys = FALSE, remove_shared_spectronaut_proteingroups = FALSE, do_plot = TRUE) {
  reset_log()
  append_log("reading Spectronaut report...", type = "info")


  rt_col = "EG.iRTEmpirical"
  if(!use_irt) rt_col = c("EG.RTEmpirical", "EG.MeanApexRT")

  ds = import_dataset_in_long_format(filename,
                                       attributes_required = list(sample_id = "R.FileName",
                                                                  protein_id = c("PG.ProteinAccessions", "PG.ProteinId", "EG.ProteinId", "FG.ProteinId", "PG.ProteinGroups"),
                                                                  sequence_modified = c("EG.ModifiedSequence", "FG.Id"),
                                                                  charge = c("FG.Charge", "FG.Id"),
                                                                  rt = rt_col,
                                                                  isdecoy = "EG.IsDecoy",
                                                                  intensity = c("FG.MS2PeakArea", "FG.TotalPeakArea", "FG.MS2RawQuantity", "FG.Quantity"),
                                                                  cscore = "EG.Cscore",
                                                                  confidence = "EG.Qvalue"),
                                       attributes_optional = list(spectral_library = "EG.Library",
                                                                  mz = "FG.PrecMzCalibrated",
                                                                  intensity_norm = c("FG.NormalizedMS2PeakArea", "FG.NormalizedTotalPeakArea", "FG.MS2Quantity")),
                                       confidence_threshold = confidence_threshold,
                                     return_decoys = return_decoys,
                                     do_plot = do_plot)

  ds$peptides = ds$peptides %>% mutate(confidence_score = cscore)

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
#' @param filename the full file path of the input file. Either supply a filename, or a data table
#' @param x input data table as a data.frame (or tibble, or data.table). Either supply a filename, or a data table
#' @param attributes_required a list of attributes that must be present in the input file (eg; sample names, peptide IDs, intensity values, etc.)
#' @param attributes_optional a list of optional attributes that we can do without downstream, like m/z
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param remove_peptides_never_above_confidence_threshold set to true to remove all target peptides that never pass the confidence score threshold
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @param select_unique_precursor_per_modseq if multiple precursors are available for a modified sequence (eg; charge 2 and 3), should we make an effort to select the best or just return all of them ?
#' @param custom_func_manipulate_DT_onload !expert use only! (there is intricate interplay between such pre-filtering and the code in this function); provide a function that manipulates the data.table after load. example implementation in the fragpipe parser.
#'
#' @importFrom data.table as.data.table chmatch setkey setorder
#' @export
import_dataset_in_long_format = function(filename = NULL, x = NULL, attributes_required=list(), attributes_optional = list(), confidence_threshold = 0.01,
                                         remove_peptides_never_above_confidence_threshold = TRUE, select_unique_precursor_per_modseq = TRUE,
                                         custom_func_manipulate_DT_onload = NULL,
                                         return_decoys = TRUE, do_plot = TRUE) {
  ### sanity check on input parameters
  if((length(filename) != 1 && length(x) == 0) || (length(filename) == 1 && length(x) > 0)) {
    append_log(paste("'filename' OR 'x' parameter is required (but don't supply both):", filename), type = "error")
  }

  # if filename is supplied, check if it's valid
  if(length(filename) == 1) {
    if(!is.character(filename) || nchar(filename) < 3) {
      append_log(paste("'filename' parameter is not valid:", filename), type = "error")
    }

    filename = path_exists(filename, NULL, try_compressed = TRUE)
    append_log(paste("input file:", filename), type = "info")
  }

  # if data table is supplied, check if it's valid
  if(length(x) > 0 && !is.data.frame(x)) {
    append_log("'x' parameter must be of type: data.frame, tibble or data.table):", type = "error")
  }
  if(length(attributes_required) == 0 || !is.list(attributes_required)) {
    append_log("'attributes_required' parameter must be a list", type = "error")
  }
  if(!is.list(attributes_optional)) {
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

  is_file = length(filename) == 1
  if (is_file) {
    # read first line from file; read 1 line into data.frame without colnames with first row having all values, then convert to character array
    # this leverages data.table to infer what character was used to delimit columns
    headers = as.character( read_textfile_compressed(filename, as_table = T, nrow = 1, header = F, data.table = F) )
  } else {
    headers = colnames(x)
  }

  # this code is also checking for required data, so run even if data table is already in memory
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)
  # collect all column indices and their desired names
  col_indices = c(map_required, map_optional)
  col_indices_dupe = duplicated(as.integer(col_indices))
  col_indices_unique = col_indices[!col_indices_dupe] # don't use unique(), it'll drop the names

  if (is_file) {
    ### try to detect data tables with a comma as decimal separation character
    # (not uncommon when using upstream raw data processing software on Windows systems for some software * system language combinations)
    dec_char = "." # initial guess
    # select two numeric columns from the data table that are guaranteed to be present, then read the first N lines of only those columns
    df_test_numeric = read_textfile_compressed(filename, nrow = 1000, skip_empty_rows = F, as_table = T, header = F, skip = 1, select = as.numeric(map_required[c("rt", "intensity")]), stringsAsFactors = F, dec = dec_char, data.table = F)
    df_test_numeric__fail = ! sapply(df_test_numeric, class, simplify = T) %in% c("numeric", "integer")
    # if any column was parsed as non-numeric with current separation char -->> try to parse with alternative separation character downstream
    if(any(df_test_numeric__fail)) {
      append_log("The input data table does not seem to use '.' for decimal notation. Trying to parse file with ',' as decimal specification instead...", type = "warning")
      dec_char = ","
    }
    rm(df_test_numeric, df_test_numeric__fail)


    ### read table
    # only read columns of interest to speed up file parsing
    # don't parse first row as header and then skip it to avoid issues with tables that have mismatched headers (e.g. MetaMorpheus files that have extra trailing tab on header row)
    DT = read_textfile_compressed(filename, skip_empty_rows = F, as_table = T, select = as.integer(col_indices_unique), header = F, skip = 1, stringsAsFactors = F, dec = dec_char, data.table = T)
    if(is.null(DT)) {
      append_log("failed to read data table from file", type = "error")
    }
    colnames(DT) = names(col_indices_unique) # overwrite column names from file with the desired names from column specification
  } else {
    x = x[,as.integer(col_indices_unique)]
    colnames(x) = names(col_indices_unique) # overwrite column names from file with the desired names from column specification

    DT = data.table::as.data.table(x) # cast input table to data.table type
  }


  ### suppose the same column from input table is assigned to multiple attributes (eg; Spectronaut FG.Id as input for both charge and modseq)
  # j = index in 'col_indices' where we should borrow content from other column
  for(j in which(col_indices_dupe)) {
    j_source_column = names(col_indices_unique)[as.integer(col_indices_unique) == as.integer(col_indices[j])]
    DT[,names(col_indices)[j]] = DT[[j_source_column]] # use column names !
  }


  ### inject function for custom manipulation of DT. example; for some input format, repair a column. for FragPipe, the modified peptide column is empty for peptides without a modification, so there is no column with proper modseq available. have to 'manually repair', but other than that we want to maintain this generic function
  if(length(custom_func_manipulate_DT_onload) == 1 && is.function(custom_func_manipulate_DT_onload)) {
    DT = custom_func_manipulate_DT_onload(DT)
  }



  ### decoys
  if("isdecoy" %in% colnames(DT)) {
    DT[ , isdecoy := isdecoy %in% c(TRUE, 1, "decoy", "true", "True", "TRUE"), by=isdecoy]
    ## apply type casting to all rows is slower than first using data.table grouping;
    # DT$isdecoy = DT$isdecoy %in% c(TRUE, 1, "decoy", "true", "True", "TRUE")
  } else {
    DT$isdecoy = FALSE
  }


  ### confidence
  # no confidence scores in the input table, assume all are 'detected'. otherwise, test confidence score vs threshold
  if(!"confidence" %in% colnames(DT)) {
    append_log("input dataset lacks peptide confidence scores, thus assuming all peptides are 'detected'", type = "warning")
    DT[ , detect := TRUE ]
  } else {
    # define detect using a threshold on confidence scores
    DT[ , detect := is.finite(confidence) & confidence <= confidence_threshold ]
  }


  ### remove entries that lack crucial data
  # 'remove' target peptides that don't have intensity >= 1. For decoys, intensity values are not required
  DT[ , rows_remove := sample_id == "" | protein_id == "" | sequence_modified  == "" | (!isdecoy & (!is.finite(intensity) | intensity < 1))]


  ### format sample_id; strip path and whitelisted extensions from filename. using grouping, we only have to apply the regex on unique elements
  DT[ , sample_id := gsub(regex_rawfile_strip_extension(), "", sample_id, ignore.case=T), by=sample_id]


  ### format charge
  # assume the charge is either provided as a plain number, or as trailing character in a string (eg; FragmentGroup ID @ spectronaut)
  if(typeof(DT$charge) == "character") { # if input data in 'charge' column is not numeric...
    DT[ , charge := replace_na(suppressWarnings(as.integer( sub(".*\\D(\\d+)$", "\\1", charge) )), 0), by=charge]
  } else {
    DT[ , charge := replace_na(suppressWarnings(as.integer(charge)), 0), by=charge]
  }

  ### format modified sequence
  # remove leading and trailing annotations from modified sequence, such as underscores or charge states. regex; start with non-character or opening bracking (modification), and analogous for trailing
  DT[ , sequence_modified := gsub("(^[^a-zA-Z([]+)|([^]a-zA-Z)]+$)", "", sequence_modified), by=sequence_modified] # relatively slow, but probably not much to be gained
  ## test data and comparison of regular expressions  +  benchmarks
  # mock = c("ABC", "_ABC", "_(n)ABC", "[n]ABC", "ABC.", "ABC.3", "ABCn3", "aa", "aa_", "_aa", "_(b)aa(c)_", "(b)aa(c)..", "_[b]aa[c]_", "2.[b]aa[c]", "_[c(unimod:1)]aaa")
  # cbind(mock, gsub("(^[^a-zA-Z([]+)|([^a-zA-Z)\\]]+$)", "", mock, perl = TRUE), gsub("(^[^a-zA-Z([]+)|([^]a-zA-Z)]+$)", "", mock, perl = F) )
  # microbenchmark::microbenchmark(
  #   a = {DT[ , sequence_modified := gsub("(^[^a-zA-Z([]+)|([^a-zA-Z)\\]]+$)", "", sequence_modified, perl = TRUE), by=sequence_modified]},
  #   b = {DT[ , sequence_modified := gsub("(^[^a-zA-Z([]+)|([^]a-zA-Z)]+$)", "", sequence_modified, perl = F), by=sequence_modified]},
  #   c = {DT[ , sequence_modified := sub("[^a-zA-Z([]*(.*)[^]a-zA-Z)]*", "\\1", sequence_modified, perl = F), by=sequence_modified]},
  #   times = 25
  # )


  # if no plain sequence is available in input table, extract it from the modified sequence by removing everything between brackets. either [] or ()
  # note that above we validated sequence_modified is never empty. so if there the provided sequence_plain has any empty values, default back to a manipulation of modified sequences
  if (!"sequence_plain" %in% colnames(DT) || any(DT$sequence_plain == "")) {
    DT[ , sequence_plain := gsub("(\\[[^]]*\\])|(\\([^)]*\\))", "", sequence_modified), by=sequence_modified] # slow !
  }


  # peptide_id = modified sequence + charge (eg. unique precursors)
  DT[ , peptide_id := paste0(sequence_modified, "_", charge), by= .(sequence_modified, charge)]

  # there should not be any duplicate peptide_id*sample_id pairs (within the subset of non-decoy rows we intend to retain)
  if(anyDuplicated(DT[isdecoy==FALSE & rows_remove==FALSE], by=c("peptide_id", "sample_id")) != 0) {
    err_message = "error: Your dataset contains multiple rows of data for some peptide_id*sample_id combinations. There should always be at most 1 row in the input data table for any combination of modified sequence + charge. Perhaps some raw file is included multiple times (e.g. C:/data/sample1.wiff and c:/data/sample1.dia   or   C:/data1/sample.raw and C:/data2/sample.raw) ?"

    if("spectral_library" %in% colnames(DT) && n_distinct(DT$spectral_library) > 1) {
      # if there is a spectral library column, perhaps the upstream software presented results per spectral library ?
      if(anyDuplicated(DT[isdecoy==FALSE & rows_remove==FALSE], by=c("peptide_id", "sample_id", "spectral_library")) != 0) {
        # throw error if there are duplicates even within spectral_library (cannot rescue; results are ambiguous)
        append_log(err_message, type = "error")
      } else {
        # otherwise, try to rescue by selecting the 'best'. Print as warning, this is unexpected input data (upstream software should decide which peptide_id to report)
        append_log("input dataset contains ambiguous peptide_id per sample spread across spectral libraries, now selecting 1 peptide_id per sample", type = "warning")

        # note; analogous to 'select_unique_precursor_per_modseq' implementation below
        # compose temp column that includes the spectral library
        DT[ , peptide_id_speclib := paste0(sequence_modified, "_", charge, "_", spectral_library), by= .(sequence_modified, charge, spectral_library)]
        # simple summary statistics per peptide_id (per spectral library)
        x = DT[rows_remove==FALSE & isdecoy==FALSE][ , .(peptide_id=peptide_id[1], ndetect=sum(detect==TRUE), abundance=mean(intensity)), by=peptide_id_speclib]
        # sort such that the desireable sequence_modified for each peptide_id comes out on top
        data.table::setorder(x, -ndetect, -abundance, na.last = TRUE)
        # select first entry
        x = x[ , index := seq_len(.N), by = peptide_id][index == 1]
        # flag for removal of peptide_id not present in subsetted table x
        DT[rows_remove==FALSE & isdecoy==FALSE & ! data.table::`%chin%`(peptide_id_speclib, x$peptide_id_speclib), rows_remove := TRUE, by=peptide_id_speclib]
        DT[ , peptide_id_speclib := NULL] # remove temp column
      }
    } else {
      append_log(err_message, type = "error")
    }
  }


  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  # note that above, we already removed (i.e. flagged as 'rows_remove') non-decoy rows where intensity values are below 1
  DT[ , intensity := suppressWarnings(log2(intensity)) ]
  DT[!is.finite(intensity) | intensity < 0, intensity := NA ] # update only rows with invalid value; threshold intensity values at 1 (but note here is log2 transformed already)
  if("intensity_norm" %in% colnames(DT) && typeof(DT$intensity_norm) == "numeric") {
    DT[ , intensity_norm := suppressWarnings(log2(intensity_norm)) ] # update entire column
    DT[!is.finite(intensity_norm) | intensity_norm < 0, intensity_norm := NA ] # update only rows with invalid value; threshold intensity values at 1 (but note here is log2 transformed already)
  }

  ### for DIA data we can compute Cscore histograms of both target and decoy peptides (which will be shown in report downstream)
  plotlist = list()
  if(do_plot && "cscore" %in% colnames(DT)) {
    plotlist$ggplot_cscore_histograms = plot_dia_cscore_histograms(as_tibble(DT))
    # bugfix; ggplot objects can have huge memory footprints, esp. when writing these objects to RData file. So we here render the plot and store only the resulting grob
    plotlist$ggplot_cscore_histograms = reduce_ggplot_object_size(plotlist$ggplot_cscore_histograms)
  }

  if(!return_decoys) {
    DT[ , rows_remove := rows_remove | isdecoy ]
  }


  ### remove peptides never detected
  # (eg; DIA dataset, some peptide from the spectral library that is never confidently identified/quantified)
  if(remove_peptides_never_above_confidence_threshold) {
    # unique set of peptide_id that have detect==TRUE at least once
    pepid_valid = DT[rows_remove==FALSE & detect==TRUE, .(peptide_id), by=peptide_id]$peptide_id
    # non-decoys that are not on remove list yet, but their peptide_id does not match valid list -->> remove
    DT[rows_remove==FALSE & isdecoy==FALSE &  ! data.table::`%chin%`(peptide_id, pepid_valid), rows_remove := TRUE, by=peptide_id] # NOT peptide_id %in% pepid_valid
  }


  ### remove peptides attributed to multiple protein identifiers
  DT_valid = DT[rows_remove==FALSE & isdecoy==FALSE][ , .(sequence_plain), by = .(sequence_plain, protein_id)][ , .( fail = .N != 1), by = sequence_plain]
  if(any(DT_valid$fail)) {
    # array of plain sequences to remove
    pseq_invalid = DT_valid[fail==TRUE, .(sequence_plain), by=sequence_plain]$sequence_plain
    append_log(sprintf("%d unique target (plain)sequences ambiguously mapped to multiple proteins and thus removed. Examples; %s",
                       length(pseq_invalid), paste(head(pseq_invalid, 5), collapse=", ")), type = "info")
    # remove non-decoy where plain sequence matches the subset that failed the above uniqueness test
    DT[rows_remove==FALSE & isdecoy==FALSE & data.table::`%chin%`(sequence_plain, pseq_invalid), rows_remove := TRUE, by=sequence_plain] # sequence_plain %chin% pseq_invalid
  }


  ### select 'best' peptide_id for each sequence_modified, in case there are multiple data points provided due to multiple charges, spectral libraries, etc.
  if(select_unique_precursor_per_modseq) {
    # # select 'best' peptide_id for each sequence_modified, in case there are multiple data points provided due to multiple charges, spectral libraries, etc.
    # data.table::setkey(DT, peptide_id)
    # simple summary statistics per peptide_id
    x = DT[rows_remove==FALSE & isdecoy==FALSE][ , .(sequence_modified=sequence_modified[1], ndetect=sum(detect==TRUE), abundance=mean(intensity)), by=peptide_id]
    npepid_pre = nrow(x)
    # sort such that the desireable sequence_modified for each peptide_id comes out on top
    data.table::setorder(x, -ndetect, -abundance, na.last = TRUE)
    # select first entry
    x = x[ , index := seq_len(.N), by = sequence_modified][index == 1]
    npepid_post = nrow(x)
    # debug; x[sequence_modified=="KNPDSQYGELIEK"]; y[sequence_modified=="KNPDSQYGELIEK"]
    # flag for removal of peptide_id not present in subsetted table x
    DT[rows_remove==FALSE & isdecoy==FALSE & ! data.table::`%chin%`(peptide_id, x$peptide_id), rows_remove := TRUE, by=peptide_id] # NOT peptide_id %in% x$peptide_id

    if(npepid_pre != npepid_post) {
      append_log(sprintf("%d/%d precursors remain after selecting the 'best' precursor for each modified sequence", npepid_post, npepid_pre), type = "info")
    }
  }

  ### finally, subset DT (do once = no garbage collection overhead)
  result = tibble_peptides_reorder(as_tibble(DT[rows_remove == FALSE]) %>% select(-rows_remove))

  if(!any(result$isdecoy == FALSE)) {
    append_log("there must be target peptides (eg; not a decoy and intensity value > 0) in the peptides tibble", type = "error")
  }

  proteins = empty_protein_tibble(result)

  # optional: log_peptide_tibble_pep_prot_counts(tib)
  return(list(peptides=result, proteins=proteins, plots=plotlist))
}
