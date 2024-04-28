
#' Import a label-free proteomics dataset from FragPipe
#'
#' @description
#' To generate output files that we here require in MS-DAP, configure FragPipe as follows:
#' 1) assign unique Experiment identifiers in the FragPipe workflow tab.
#' Importantly, Fragpipe produces output (that MS-DAP can import) for each unique Experiment*Bioreplicate combination.
#' So if you provide the same experiment ID for all input files, the output looks like 1 sample to MS-DAP which is not what we want;
#' ensure all experiments have unique values in the FragPipe workflow tab !
#' 2) enable IonQuant in the "Quant (MS1)" tab
#' (typically, this is already done if you loaded the LFQ-MBR workflow in the FragPipe workflow tab).
#' 3) optionally, enable match-between-runs in the FragPipe "Quant (MS1)" tab.
#'
#' This function will merge data from various FragPipe output files to construct a MS-DAP dataset.
#' If you followed above instructions, FragPipe should have generated the following files in the
#' output directory that will be parsed by this function;
#' - combined_protein.tsv
#' - for each experiment, a folder that includes 2 MS-DAP required files: ion.tsv and psm.tsv
#'
#' (as most recently tested with FragPipe v18 and 19.1)
#'
#' @param path the full file path to the FragPipe output directory
#' @param acquisition_mode the type of experiment, should be a string. Valid options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', should be a numeric value between 0 and 1 (peptides with a 'Expectation' value in psm.tsv that is <= this threshold are classified as 'detected')
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_fragpipe_ionquant = function(path, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()

  ### check for valid function parameters
  check_parameter_is_string(path)
  if(!dir.exists(path)) {
    append_log(paste("directory does not exist:", path), type = "error")
  }

  check_parameter_is_string(acquisition_mode)
  if( ! acquisition_mode %in% c("dda","dia")) {
    append_log('acquisition_mode parameter must be either "dda" or "dia"', type = "error")
  }

  check_parameter_is_numeric(confidence_threshold, minval_ = 0, maxval_ = 1)
  check_parameter_is_string(collapse_peptide_by)
  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }



  ### check whether all required input files are present
  append_log("reading FragPipe output folder", type = "info")
  errmsg = 'Did you provide experiment information in the FragPipe Workflow tab? If you do not enter experiment info (either unique combinations of "Experiment" and "Bioreplicate", or unique values for "Experiment" throughout), FragPipe does not generate required output files for MS-DAP. See further online documentation at MS-DAP GitHub > user guide > FragPipe'

  file_combined_protein = path_exists(path, "combined_protein.tsv", try_compressed = TRUE, silent = TRUE)
  if(file_combined_protein == "") {
    append_log(sprintf('Cannot find file "combined_peptide.tsv" at provided path "%s"\n%s', path, errmsg), type = "error")
  }

  subdirs = list.dirs(path, full.names = TRUE, recursive = FALSE)
  if(length(subdirs) == 0) {
    append_log(sprintf('Cannot find any subdirectories, which should contain the ion.tsv and psm.tsv FragPipe result files, at provided path "%s"\n%s', path, errmsg), type = "error")
  }

  files_psm = files_ion = NULL
  for(d in subdirs) {
    fpsm = path_exists(d, "psm.tsv", try_compressed = TRUE, silent = TRUE)
    fion = path_exists(d, "ion.tsv", try_compressed = TRUE, silent = TRUE)
    if(fpsm == "") {
      append_log(sprintf('Cannot find "psm.tsv" file in FragPipe results at "%s"\n%s', d, errmsg), type = "error")
    }
    if(fion == "") {
      append_log(sprintf('Cannot find "ion.tsv" file in FragPipe results at "%s"\n%s', d, errmsg), type = "error")
    }
    files_psm = c(files_psm, fpsm)
    files_ion = c(files_ion, fion)
  }

  ## optional file that we can use to double-check FragPipe configuration
  ## note; check disable as it's redundant with downstream checks in psm.tsv
  # file_experiment_annotation = path_exists(path, "experiment_annotation.tsv", try_compressed = TRUE, silent = TRUE)
  # if(file_experiment_annotation != "") {
  #   # read tsv table
  #   x = as_tibble(read_textfile_compressed(file_experiment_annotation, as_table = T))
  #   # if it contains columns condition and replicate, we can verify that the combination thereof is unique
  #   if(all(c("condition", "replicate") %in% colnames(x))) {
  #     id = paste(x$condition, x$replicate)
  #     if(anyDuplicated(id)) {
  #       append_log("'experiment_annotation.tsv' describes multiple samples were assigned the same experiment*replicate value; this should only occur when your dataset was fractionated (thus multiple raw files per 'sample'). If your dataset does not contain fractionated samples, this indicates FragPipe was not configured correctly for use with MS-DAP (in the FragPipe workflow tab, ensure each sample has a unique combination for Experiment*Bioreplicate)", type = "warning")
  #     }
  #   }
  #   rm(x)
  # }



  ### iterate all psm and ion files, parsing respective file pairs
  tib_result = NULL
  for(i in seq_along(files_psm)) {
    tib_result = bind_rows(
      tib_result,
      parse_fragpipe_psm_ion_pair(file_psm = files_psm[i], file_ion = files_ion[i])
    )
  }


  ### updated protein_id with ambiguous included
  protid_to_fullprotid = parse_fragpipe_proteins_from_combined_proteins(file_combined_protein)
  # ensure the protein table covers all proteins reported in psm.tsv
  pid_mismatch = setdiff(unique(tib_result$protein_id), protid_to_fullprotid$protein_id)
  if(length(pid_mismatch) > 0) {
    append_log(sprintf('%d unique protein ID from psm.tsv and ion.tsv files could not be matched against combined_protein.tsv , so "Indistinguishable Proteins" are not available for these protein groups (does not affect MS-DAP data processing, only affects assigned protein and gene ID for respective protein groups which now only describe the leading protein)', length(pid_mismatch)), type = "warning")
  }
  # update protein_id in peptide table
  tib_result = tib_result %>%
    left_join(protid_to_fullprotid, by = "protein_id") %>%
    # if we found the protein_id in the master table, update protein-group-ID
    mutate(protein_id = ifelse(is.na(protein_id_full), protein_id, protein_id_full)) %>%
    select(-protein_id_full)


  ### format peptide table columns, e.g. ensure intenties are log-transformed
  tib_result = tib_result %>%
    mutate(
      detect = is.finite(confidence) & confidence <= confidence_threshold, # define 'detect'
      rt = rt / 60, # convert RT to minutes
      intensity = log2(intensity), # log2 transform intensities
    )


  ### collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  log_peptide_tibble_pep_prot_counts(tib_result)
  # report MBR hits based on 'has RT'. Although we could use the confidence value / detect flag, this allows us to monitor if we pulled any values from ion.tsv that are not in psm.tsv
  tmp = sum(!is.finite(tib_result$rt))
  append_log(sprintf("%d peptide_id * sample_id in the peptide table imported from FragPipe (after postprocessing), of which %d (%.1f%%) are MBR hits", nrow(tib_result), tmp, tmp/nrow(tib_result)*100), type = "info")
  return(list(peptides=tibble_peptides_reorder(tib_result), proteins=empty_protein_tibble(tib_result), plots=list(), acquisition_mode = acquisition_mode))
}



#' augment FragPipe ion.tsv with retention-time data from the respective psm.tsv file
#'
#' note that this still leaves retention-time information for MBR hits, this data is not available in FragPipe output atm.
#'
#' @param file_ion full path to a ion.tsv file
#' @param file_psm full path to a psm.tsv file, this must match the `file_ion` argument (e.g. files in the same dir)
parse_fragpipe_psm_ion_pair = function(file_ion, file_psm) {
  spectrum = NULL
  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  append_log(paste("reading FragPipe file;", file_ion), type = "info")
  ion = read_table_by_header_spec(
    file = file_ion,
    attributes_required = list(
      protein_id = "Protein",
      sequence_plain = "Peptide Sequence",
      sequence_modified = "Modified Sequence",
      charge = "Charge",
      mz_theoretical = "M/Z",
      intensity = "Intensity",
      confidence = "Expectation"
    ),
    as_tibble_type = T
  ) %>%
    # remove ions without intensity values
    filter(is.finite(intensity) & intensity > 0) %>%
    mutate(
      protein_id = fasta_id_short(protein_id),
      peptide_id = paste(sequence_modified, charge, sep = "_"),
      confidence = ifelse(is.finite(confidence), confidence, NA)
    ) %>%
    # don't assume the modified sequence + charge are unique
    # (should be, but deal with it anyway in case FragPipe updates something upstream)
    group_by(peptide_id) %>%
    summarise(
      protein_id = protein_id[1],
      sequence_plain = sequence_plain[1],
      sequence_modified = sequence_modified[1],
      charge = charge[1],
      mz_theoretical = mz_theoretical[1],
      intensity = sum(intensity),
      confidence = suppressWarnings(min(confidence, na.rm = T)), # can be all missing / NA (e.g. MBR)
      .groups = "drop"
    ) %>%
    mutate(
      confidence = ifelse(is.finite(confidence), confidence, NA),
      # temp vars that we'll repeatedly use downstream
      is_mbr = !is.finite(confidence),
      is_modified = sequence_plain != sequence_modified,
      key_plainseq_charge = paste(sequence_plain, charge)
    )

  if(nrow(ion) == 0) {
    append_log(paste("empty table;", file_ion), type = "error")
  }

  append_log(paste("reading FragPipe file;", file_psm), type = "info")
  psm = read_table_by_header_spec(
    file = file_psm,
    attributes_required = list(
      spectrum = "Spectrum",
      filename = "Spectrum File",
      sequence_plain = "Peptide",
      sequence_modified = "Modified Peptide",
      charge = "Charge",
      mz_theoretical = "Calculated M/Z",
      confidence = "Expectation",
      rt = "Retention"
    ),
    as_tibble_type = T
  ) %>%
    # temp vars that we'll repeatedly use downstream
    mutate(
      is_modified = sequence_modified != "",
      key_plainseq_charge = paste(sequence_plain, charge)
    )

  if(nrow(psm) == 0) {
    append_log(paste("empty table;", file_psm), type = "error")
  }


  # specifically for fragpipe; strip full path + "interact-" + ".pep.xml" extension
  # example entry in input table; C:\DATA\fragpipe_test\interact-a05191.pep.xml
  ufilename = sub("^interact[^a-z](.*)\\.pep\\.xml$", "\\1", basename(unique(psm$filename)))
  # remove filename column from PSM table
  psm = psm %>% select(-filename)
  # sample_id for this set of raw files
  combined_sample_id = paste(unique(ufilename), collapse=";") # take unique again after removing path/basename
  # if there are many files (e.g. sample with 25 fractions and long filenames) and thus we produce a really long sample_id, create a trimmed variant
  if(length(ufilename) > 1 && nchar(combined_sample_id) > 200) {
    combined_sample_id = paste0(ufilename[1], ";", length(ufilename)-1, "-other-files")
  }

  if(length(ufilename) != 1) {
    append_log(sprintf("multiple (%d) raw files are described in a single psm.tsv file. This should only occur when your dataset was fractionated (thus multiple raw files per 'sample'). If your dataset does not contain fractionated samples, this indicates FragPipe was not configured correctly for use with MS-DAP (in the FragPipe workflow tab, ensure each sample/rawfile has a unique combination for Experiment*Bioreplicate)", length(ufilename)), type = "warning")
  }


  ### fetch peptide retention times from PSM file
  #
  ## matching between psm.tsv and ion.tsv is not straight-forward due to 2 FragPipe issues;
  # 1) cannot match directly by modified sequence, different values for modified are used in psm.tsv and ion.tsv
  # example; psm.tsv = n[43]AADTQVSETLK    ion.tsv = n[42.0106]AADTQVSETLK
  # 2) there are rounding issues in the FragPipe data, values between psm.tsv and ion.tsv sometimes don't match
  # example; Calculated Peptide Mass; psm.tsv = 1239.632   ion.tsv = 1239.6329
  #
  ## workaround for finding the closest match. Our strategy is chosen to
  ## A) not make any assumptions on the nature of these numeric mismatches (e.g. precision, or rounding, and how to round)
  ## B) avoid computations that are generally slow in R (e.g. looping over all rows in the ion table and matching to psm table with some tolerance)
  # pseudocode for matching modified sequences;
  # 1. join using a combination of plain-sequence and charge (may lead to 1:n matches)
  # 2. compute the difference in numeric value (e.g. peptide mass) between ion.tsv and psm.tsv
  # 3. sort the table and take closest match
  # (we can vectorize all this, so it only costs a few extra joins, matches and some RAM)
  #
  ## note that this still leaves all MBR hits without a retention time, which is unavoidable as this information is currently missing in FragPipe output


  # first, match non-modified sequences
  ion_nomod_nombr = ion %>%
    select(peptide_id, key_plainseq_charge, is_modified, is_mbr) %>%
    filter(is_modified == FALSE & is_mbr == FALSE) %>%
    left_join(
      # from the PSM table, take the RT of the highest confidence score PSM (of each 'plain-sequence * charge' combination)
      psm %>%
        filter(is_modified == FALSE) %>%
        arrange(confidence) %>%
        distinct(key_plainseq_charge, .keep_all = TRUE) %>%
        select(key_plainseq_charge, rt),
      by = "key_plainseq_charge") %>%
    select(peptide_id, rt) %>%
    filter(is.finite(rt)) # only retain peptide_id where our RT lookup succeeded


  ## find closest match per modified sequence
  # importantly, we join by plainsequence*charge first thus ordering by mz_diff finds the closest precursor in terms of "theoretical m/z provided by FragPipe" from only peptides with the exact same plain sequence and charge
  # 1) the 'Observed Mass' columns in psm.tsv and ion.tsv don't align, so we use the theoretical M/Z instead
  # 2) because we select from ion.tsv only modified sequences that are not match-between-runs, we expect each precursor to be present in psm.tsv
  ion_mod_nombr = ion %>%
    select(peptide_id, key_plainseq_charge, mz_theoretical_ion = mz_theoretical, is_modified, is_mbr) %>%
    filter(is_modified == TRUE & is_mbr == FALSE) %>%
    # join by plain sequence*charge, so 1:n join is possible (intentionally, we don't know which modified sequence to match due to FragPipe unstable sequence IDs)
    left_join(
      psm %>%
        filter(is_modified == TRUE) %>%
        select(spectrum, key_plainseq_charge, mz_theoretical_psm = mz_theoretical, confidence, rt),
      by = "key_plainseq_charge"
    ) %>%
    filter(is.finite(rt)) %>% # only retain peptide_id where our RT lookup succeeded
    mutate(mz_diff = abs(mz_theoretical_psm - mz_theoretical_ion)) %>%
    # take closest match, if multiple (i.e. sorting ties) then take highest-confidence PSM
    arrange(mz_diff, confidence) %>%
    select(peptide_id, rt, mz_diff) %>%
    distinct(peptide_id, .keep_all = TRUE)


  ### map results back into ion table
  ion = ion %>%
    select(peptide_id, protein_id, sequence_plain, sequence_modified, charge, confidence, intensity) %>%
    mutate(
      # Note that we simply set the same sample filename for all PSM (first filename encountered in psm.tsv) because
      # the ion.tsv file does not contain information to trace back to the raw file where the peakarea was derived from.
      # So the peak area ion ion.tsv could be the sum over multiple files, or a MBR hit not listed in psm.tsv at all.
      # Possible solutions are to just assign 1 filename, or paste all listed filenames together as 1 'id'
      sample_id = combined_sample_id,
      # init column
      rt = NA_real_
    )

  if(nrow(ion_nomod_nombr) > 0) {
    ion$rt[data.table::chmatch(ion_nomod_nombr$peptide_id, ion$peptide_id)] = ion_nomod_nombr$rt
  }

  if(nrow(ion_mod_nombr) > 0) {
    ion$rt[data.table::chmatch(ion_mod_nombr$peptide_id, ion$peptide_id)] = ion_mod_nombr$rt
    # check if any modified peptides were matched by large m/z diff and issue warning
    tmp = tail(ion_mod_nombr$mz_diff, 1) # already sorted by m/z diff
    if(tmp > 0.001) {
      append_log(sprintf("searching for each precursor's retention-time; matching between FragPipe ion.tsv and psm.tsv, largest m/z mismatch within peptides of the same (plain) sequence was %s", tmp), type = "info")
    }
  }

  return(ion)
}



#' Construct MS-DAP style protein group identifiers from a FragPipe 'combined_protein.tsv' file by appending all 'indistinguishable' proteins after the leading protein ID
#'
#' @param filename full path to a FragPipe 'combined_protein.tsv' file
parse_fragpipe_proteins_from_combined_proteins = function(filename) {
  x = as_tibble(read_textfile_compressed(filename, as_table = T))
  # input validation
  col_missing = setdiff(c("Protein", "Indistinguishable Proteins"), colnames(x))
  if (length(col_missing) > 0) {
    append_log(sprintf('Expected column names in %s table: "Protein", "Indistinguishable Proteins". Missing: %s', filename, paste(col_missing, collapse = ", ")), type = "error")
  }

  x = x %>%
    select(protein_id = Protein, protein_id_ambiguous = `Indistinguishable Proteins`) %>%
    # enforced unique `protein_id`
    distinct(protein_id, .keep_all = TRUE) %>% # shouldn't be needed, but enforce unique protein IDs in the table we just read from file
    mutate(
      protein_id = fasta_id_short(protein_id), # enforce short protein IDs
      protein_id_ambiguous_list = strsplit(protein_id_ambiguous, "[ ,;]+") # to list column
    )

  y = x %>%
    # to long format
    unnest(cols = protein_id_ambiguous_list, keep_empty = FALSE) %>%
    # now we can apply vectorized function to enforce short protein IDs
    mutate(protein_id_ambiguous_list = fasta_id_short(protein_id_ambiguous_list)) %>%
    # back to protein-group format
    group_by(protein_id) %>% # enforced unique `protein_id`
    summarise(protein_id_ambiguous = paste(protein_id_ambiguous_list, collapse = ";"), .groups = "drop")

  x %>%
    select(protein_id) %>%
    left_join(y, by = "protein_id") %>%
    mutate(protein_id_full = ifelse(is.na(protein_id_ambiguous), protein_id, paste(protein_id, protein_id_ambiguous, sep = ";"))) %>%
    select(-protein_id_ambiguous)
}
