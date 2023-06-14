
#' Import a label-free DDA proteomics dataset from a ProteomeDiscoverer PSM result file
#'
#' @description
#' ProteomeDiscoverer workflow must include Percolator so MS-DAP can parse peptide confidence scores.
#'
#' Example ProteomeDiscoverer workflow:
#' - Processing Step: PWF_QE_Precursor_Quan_and_LFQ_SequestHT_Percolator
#' - Consensus Step: CWF_Comprehensive_Enhanced Annotation_LFQ_and_Precursor_Quan
#' - Consensus Step: add the "result exporter" (drag&drop from side panel to bottom panel)
#'
#'
#' Optionally, you can relax the output filter criteria somewhat and use MS-DAP filtering as follows;
#' - Consensus step -->> "peptide and protein filter" -->> Peptide Confidence At Least -->> change to medium
#' - While you're making changes there, you can also set "Remove Peptides Without Proteins" to "True"
#'
#'
#' Besides validating the input data table and reformatting table columns for MS-DAP, this function;
#' 1) for each precursor (modified peptide sequence + charge state), the best PSM confidence score is retained.
#' 2) modified peptide sequences are reformatted from input table where sequence and modifications are stored separately.
#' e.g. Annotated Sequence = "AAGLATmISTmRPDIDNmDEYVR" and Modifications = "M7(Oxidation); M11(Oxidation); M18(Oxidation)" in input becomes modified_sequence in MS-DAP; "AAGLATM(Oxidation)ISTM(Oxidation)RPDIDNM(Oxidation)DEYVR"
#'
#' @param filename full path to the ProteomeDiscoverer _PSMs.txt file
#' @param confidence_threshold confidence score threshold ('Percolator q-Value' column in the PSM file) at which a peptide is considered 'identified', default: 0.01 (target value must be lesser than or equals)
#' @param remove_lowconf boolean value indicating whether peptides classified as 'low confidence' by ProteomeDiscoverer should be removed from the results
#' @param one_psm_per_precursor optionally, retain for each precursor in each sample only the peakarea for 1 PSM.
#' This parameter allows you to control how abundance values from precursors matched by multiple PSM are handled, as this might depend on your ProteomeDiscoverer settings.
#' If ProteomeDiscoverer performed peak integration and reports the same (redundant) peak intensity for each PSM of the same precursor, we suggest to use `one_psm_per_precursor = "intensity"`.
#' Note that relevant statitics for your dataset will be printed to the console/log (e.g. fraction of redundant PSM that contain unique intensity values).
#' Set to "" to disable this filtering and use the sum of all PSM intensity values per precursor*sample (default).
#' Use `one_psm_per_precursor = "intensity"` to select the highest intensity value (within the subset of PSM where `confidence < confidence_threshold`).
#' Use `one_psm_per_precursor = "confidence"` to select the intensity value from the PSM with best/lowest confidence value
#' @param collapse_peptide_by if multiple data points are available for a peptide sequence in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_proteomediscoverer_txt = function(filename, confidence_threshold = 0.01, remove_lowconf = TRUE, one_psm_per_precursor = "", collapse_peptide_by = "sequence_modified") {
  reset_log()
  append_log("reading ProteomeDiscoverer PSM file", type = "info")

  check_parameter_is_string(filename)
  check_parameter_is_numeric(confidence_threshold)
  check_parameter_is_boolean(remove_lowconf)
  check_parameter_is_string(one_psm_per_precursor)
  check_parameter_is_string(collapse_peptide_by)

  if(!(one_psm_per_precursor %in% c("confidence", "intensity", ""))) {
    append_log('one_psm_per_precursor parameter must be any of; "confidence", "intensity", ""', type = "error")
  }
  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }


  # will check for presence of file as well as .gz/.zip extension if file doesn't exist, will throw error if both do not exist
  filename = path_exists(filename, try_compressed = TRUE)

  ### read evidence file and validate input columns
  attributes_required = list(sequence_fullstring = "Annotated Sequence",
                             sequence_modificationstring = "Modifications",
                             protein_id = "Master Protein Accessions",
                             intensity = "Precursor Abundance",
                             proteingroup_count = c("# Protein Groups", "Number of Protein Groups"),
                             rt = c("Apex RT [min]", "Apex RT in min", "RT [min]", "RT in min"),
                             raw_file = "Spectrum File",
                             confidence = c("Percolator q-Value"), # use the percolator q-value as confidence metric
                             charge = "Charge")
  attributes_optional = list(mz = c("m/z [Da]", "mz in Da"),
                             confidence_classification = "Confidence")

  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  x = read_table_by_header_spec(filename, attributes_required, attributes_optional, as_tibble_type = TRUE)

  # force numeric types
  x$confidence = suppressWarnings(as.numeric(x$confidence))
  x$proteingroup_count = suppressWarnings(as.integer(x$proteingroup_count))
  x$rt = suppressWarnings(as.numeric(x$rt))
  x$charge = suppressWarnings(as.integer(x$charge))
  x$intensity = suppressWarnings(as.numeric(x$intensity))


  ### remove erronous rows;
  rows_fail = !is.finite(x$confidence) | !is.finite(x$proteingroup_count) | !is.finite(x$rt) | !is.finite(x$charge) | !is.finite(x$intensity) |
    x$sequence_fullstring == "" | x$protein_id == "" | x$raw_file == "" | x$intensity <= 0
  nrow_missingdata = sum(rows_fail)
  rows_fail = rows_fail | x$proteingroup_count != 1L
  nrow_ambiguous = sum(rows_fail) - nrow_missingdata
  nrow_total = nrow(x)
  append_log(sprintf("%d/%d rows in PSM file will be removed because of missing data (e.g. PSM without intensity value), %d/%d rows will be removed because the peptide was assigned to multiple proteingroups (PSM table from PD, column '# Protein Groups')",
                     nrow_missingdata, nrow_total, nrow_ambiguous, nrow_total), type = "info")

  # optionally, remove peptides classified as low-confidence by PD @ PSM table column 'confidence'
  if(remove_lowconf) {
    if(!"confidence_classification" %in% colnames(x)) {
      append_log("parameter 'remove_lowconf' specified the removal of peptides classified as low-confident by ProteomeDiscoverer, but the 'Confidence' column was missing from the input PD PSM table", type = "error")
    }
    rows_fail = rows_fail | tolower(x$confidence_classification) == "low"
  }

  # finally, subset data table
  x = x[rows_fail == FALSE, ]


  ### QC
  # peptide intensity = sum of peak intensities reported for respective PSM
  # This is computed downstream @ peptides_collapse_by_sequence()
  # However, it seems that in some (but not all!) PD configurations, peak integration is performed at precursor level
  # and thus the input table for this function contains the same "intensity" value for all PSM of the same precursor.
  # In this case, one should not sum abundance values as that will cause bias / intensity inflation for precursors with multiple PSM.
  # We here quantify this and log relevant statistics so users can make informed choices
  # id that combines precursor*sample
  x$groupid_temp = paste(x$sequence_fullstring, x$sequence_modificationstring, x$charge, x$raw_file)
  groupid_temp_multi = unique(x$groupid_temp[duplicated(x$groupid_temp)])

  if(length(groupid_temp_multi) > 10) {
    # for efficiency, don't check entire dataset
    if(length(groupid_temp_multi) > 10000) {
      set.seed(123)
      groupid_temp_multi = sample(groupid_temp_multi, size = 10000, replace = FALSE) # grab 5000 random entries
    }

    # subset of data table x for N unique IDs that are potentially present with redundant intensity values
    x_subset = x %>% select(groupid_temp, intensity) %>% filter(groupid_temp %in% groupid_temp_multi)
    # for each precursor*sample with multiple PSM, test if there are any unique intensity values (absolute differences as compared to the first value is larger than X)
    x_subset_test = x_subset %>% group_by(groupid_temp) %>% summarise(has_unique_values = any(abs(intensity[-1] - intensity[1]) > 0.01), .groups = "drop")

    precursor_test_count = nrow(x_subset_test)
    precursor_test_ratio_unique = sum(x_subset_test$has_unique_values) / nrow(x_subset_test)

    append_log(sprintf("while checking %d precursors were ProteomeDiscoverer reports multiple PSM, %d/%d (%.1f%%) of precursors contain unique intensity values across PSM (within the same precursor). If this is low (e.g. < 10%%), ProteomeDiscoverer is reporting the same intensity value for all PSM per precursor and you should be using the parameter 'one_psm_per_precursor==\"intensity\"' while importing your ProteomeDiscoverer dataset using MS-DAP",
                       nrow(x_subset_test), sum(x_subset_test$has_unique_values), nrow(x_subset_test), sum(x_subset_test$has_unique_values) / nrow(x_subset_test) * 100), type = "info")

    x$groupid_temp = NULL
    rm(groupid_temp_multi, x_subset, x_subset_test)
  }


  if(one_psm_per_precursor == "" && precursor_test_ratio_unique < 0.1 && precursor_test_count > 100) {
    append_log("precursors identified by multiple PSM seem to contain the same intensity value per PSM (i.e. ProteomeDiscoverer already integrated peakareas across PSM for the same precursor). Consider changing this function's parameter 'one_psm_per_precursor' from \"\" to \"intensity\" or \"confidence\"", type = "warning")
  }


  ### filter 1 PSM per precursor
  # optionally, per peptide_id*sample_id we only retain the PSM with highest confidence score
  if(one_psm_per_precursor == "confidence") {
    # note; 'sequence_fullstring' by itself are not unique (e.g. the actual modification per AA is ambiguous)
    x$peptide_id_temp = paste(x$sequence_fullstring, x$sequence_modificationstring, x$charge)
    # retain 1 row per peptide_id per sample
    x = x %>%
      arrange(confidence) %>%
      distinct(peptide_id_temp, raw_file, .keep_all = T) %>%
      select(-peptide_id_temp)
  }


  # optionally, per peptide_id*sample_id we only retain the PSM with highest abundance score (among subset that meets confidence cutoff)
  # analogous to above
  if(one_psm_per_precursor == "intensity") {
    x$peptide_id_temp = paste(x$sequence_fullstring, x$sequence_modificationstring, x$charge)
    x = x %>%
      arrange(desc(confidence <= confidence_threshold), desc(intensity)) %>%
      distinct(peptide_id_temp, raw_file, .keep_all = T) %>%
      select(-peptide_id_temp)
  }



  ### check for the input error / PD misconfiguration that yields 'ProteinCenter:sp_canonical' as protein ID instead of actual uniprot accessions
  # instead of checking for this known 'error string', we just count if there are any protein identifiers that are present suspeciously often
  uprot_qc = sort(table(x$protein_id), decreasing = T)
  uprot_qc__extreme = uprot_qc / nrow(x) # the fraction of rows in the input table that contains each protein_id
  uprot_qc__extreme_index = which(uprot_qc__extreme >= 0.1)
  if(length(uprot_qc__extreme_index) > 0) {
    append_log(paste0("at least 25% of ProteomeDiscoverer PSMs are assigned to these protein identifiers;\n", paste(names(uprot_qc__extreme)[uprot_qc__extreme_index], collapse = "\n")), type = "warning")
  }


  ### update filenames; remove extensions
  # strip path and whitelisted extensions from filename
  fnames = x %>% distinct(raw_file) %>% mutate(sample_id = gsub(regex_rawfile_strip_extension(), "", raw_file, ignore.case=T))
  x$sample_id = fnames$sample_id[data.table::chmatch(x$raw_file, fnames$raw_file)]


  ### collect all unique (modified) sequences
  # [K].AAGLATmISTmRPDIDNmDEYVR.[N]	M7(Oxidation); M11(Oxidation); M18(Oxidation)
  # [K].AAGLATMISTmRPDIDNmDEYVR.[N]	M11(Oxidation); M18(Oxidation)
  seq = x %>% select(sequence_fullstring, sequence_modificationstring) %>% distinct(sequence_fullstring, .keep_all = T) %>%
    mutate(
      # strip the bracket notation for prev/next AA
      sequence_modified = gsub("^.*\\]\\.|\\.\\[.*$|\\.$|^\\.", "", sequence_fullstring),
      # at this point, the modified AA are denoted only by lower-case but the exact modification is not included in the sequence
      sequence_plain = toupper(sequence_modified)
    )


  ### reformat modified sequence notations
  # example table row; sequence_modified = "YmHSGPVVAmVWEGLNVVK", sequence_modificationstring = "M2(Oxidation); M10(Oxidation)"
  modseq = seq %>% filter(sequence_modificationstring != "")
  if(nrow(modseq) > 0) {
    l = strsplit(modseq$sequence_modificationstring, " *; *")
    ul = data.frame(row = rep(seq_along(l), lengths(l)), input = unlist(l)) # row = index in 'modseq' table
    ul$newaastring = sub("^([a-zA-Z]+)(\\d+)(.*)", "\\1\\3", ul$input)
    ul$col = suppressWarnings(as.integer(sub("^([a-zA-Z]+)(\\d+)(.*)", "\\2", ul$input))) # col = character index within the sequence string
    # if there is no AA index in the expected/standardized modification notation, assume it's a C/N-terminal modification
    ul$fail = !is.finite(ul$col)
    if(any(ul$fail)) {
      # PD typical notation for an N-terminal modification is; N-Term(Prot)(Met-loss+Acetyl)
      ul$isnterm = ul$iscterm = F
      ul$isnterm[ul$fail] = grepl("N([ -]){0,1}Term", ul$newaastring[ul$fail], ignore.case = T)
      ul$iscterm[ul$fail] = grepl("C([ -]){0,1}Term", ul$newaastring[ul$fail], ignore.case = T)
      # default; prefix modification string provided by PD as-is  (fallback for modifications not anticipated in this codebase)
      # replace first AA with the AA prefixed with the modification
      rows = ul$fail & (ul$isnterm | !ul$iscterm)
      if(any(rows)) {
        ul$col[rows] = 1L
        ul$newaastring[rows] = paste0(gsub("[CN]([ -]){0,1}Term|\\({0,1}Prot\\){0,1}", "", ul$newaastring[rows], ignore.case = T, perl = T),
                                      substring(modseq$sequence_plain[ul$row[rows]], 1, 1) ) # first character in respective sequence
      }
      # replace last AA with the AA appended with the modification
      rows = ul$fail & ul$iscterm
      if(any(rows)) {
        ul$col[rows] = nchar(ul$newaastring[rows])
        ul$newaastring[rows] = paste0(stringr::str_sub(modseq$sequence_plain[ul$row[rows]], -1),  # last character in respective sequence
                                      gsub("[CN]([ -]){0,1}Term|\\({0,1}Prot\\){0,1}", "", ul$newaastring[rows], ignore.case = T))
      }
    }

    # convert modified sequences to a matrix where each column is an AA
    aamatrix = stringr::str_split_fixed(modseq$sequence_modified, pattern = "", n = Inf)
    # vectorized replacement; map the row/col to indices in aamatrix
    ul$matindex = ul$row + (ul$col - 1L) * nrow(aamatrix)

    ## collapse the updated modified sequences in the character matrix
    aamatrix[ul$matindex] = ul$newaastring
    modseq$sequence_modified = apply(aamatrix, 1, paste, collapse="")

    ## map these back to the input table 'seq' that holds all sequences
    i = data.table::chmatch(seq$sequence_fullstring, modseq$sequence_fullstring) # match using the raw/unformatted sequences
    seq$sequence_modified[!is.na(i)] = modseq$sequence_modified[na.omit(i)]
    # QC; print input and updated sequences;  print(modseq, n = 100)
  }

  ## the 'seq' table now contains cleaned peptide sequences; map these back to the complete peptide input tibble
  i = data.table::chmatch(x$sequence_fullstring, seq$sequence_fullstring) # match using the raw/unformatted sequences
  stopifnot(!is.na(i))
  x$sequence_plain = seq$sequence_plain[i]
  x$sequence_modified = seq$sequence_modified[i]


  ### finally, format table for downstream compatability
  x$peptide_id = paste0(x$sequence_modified, "_", x$charge)
  x$isdecoy = F
  x$intensity = log2(x$intensity)
  x$intensity[!is.na(x$intensity) & x$intensity < 0] = 0 # note; we already removed zero intensity values when importing. here, we threshold extremely low values
  x$detect = x$confidence <= confidence_threshold # peptide-level FDR cutoff

  # retain selected columns
  x = x %>% select(peptide_id, protein_id, sample_id, sequence_plain, sequence_modified, confidence, detect, rt, intensity, any_of(c("mz", "confidence_classification")))

  # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    x = peptides_collapse_by_sequence(x, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    x = peptides_collapse_by_sequence(x, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  log_peptide_tibble_pep_prot_counts(x)
  return(list(peptides = tibble_peptides_reorder(x), proteins = empty_protein_tibble(x), acquisition_mode = "dda"))
}
