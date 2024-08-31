
#' Import a label-free proteomics dataset from FragPipe; combines quantitative data from the MSstats.csv file with PSM data from psm.tsv files
#'
#' @description
#' This function requires the following FragPipe output data:
#'
#' - 'combined_protein.tsv' file, located in the FragPipe output folder
#' - 'MSstats.csv' file, located in the FragPipe output folder
#' - 'psm.psv' files, located in subdirectories of the FragPipe output folder (these are datamined to obtain peptide PSM confidence and retention times)
#'
#' Intensity values for each precursor (modified sequence and charge, columns "PeptideSequence" and "PrecursorCharge")
#' in each sample ("Run" column) are extracted from the MSstats.csv file. Next, the psm.tsv files are parsed to obtain
#' retention times and PSM confidence values for each precursor*sample.
#'
#' Unfortunately, retention times at apex peak aren't readily available for FragPipe Ionquant results across FragPipe
#' versions so in this function we obtain peptide retention times from PSM matches for now. While this yield less
#' accurate RT values and misses RT values for MBR hits (i.e. there is no PSM), the data required for this approach
#' is available for all FragPipe versions since at least 2020.
#'
#' Finally, the combined_protein.tsv file is used to obtain ambiguous protein IDs per proteingroup
#' (column "Indistinguishable Proteins").
#'
#' @param path the full file path to the FragPipe output directory
#' @param acquisition_mode the type of experiment, should be a string. Valid options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified', should be a numeric value between 0 and 1 (peptides with a '1 - PeptideProphet.Probability' value in psm.tsv that is <= this threshold are classified as 'detected')
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_fragpipe_ionquant = function(path, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()

  #### check for valid function parameters
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


  #### check whether all required input files are present
  append_log("reading FragPipe output folder", type = "info")
  errmsg = 'Did you provide experiment information in the FragPipe Workflow tab? If you do not enter experiment info (either unique combinations of "Experiment" and "Bioreplicate", or unique values for "Experiment" throughout), FragPipe does not generate required output files for MS-DAP. See further online documentation at MS-DAP GitHub > user guide > FragPipe'

  file_msstats = path_exists(path, "MSstats.csv", try_compressed = TRUE, silent = TRUE)
  if(file_msstats == "") {
    append_log(sprintf('cannot find file "MSstats.csv" at provided path "%s"\n%s', path, errmsg), type = "error")
  }

  file_combined_protein = path_exists(path, "combined_protein.tsv", try_compressed = TRUE, silent = TRUE)
  if(file_combined_protein == "") {
    append_log(sprintf('cannot find file "combined_protein.tsv" at provided path "%s"\n%s', path, errmsg), type = "error")
  }


  #### parse MSstats.csv table
  quant = fragpipe_parse_msstats(file_msstats)
  msstats_sampleid_unique = unique(quant$sample_id)
  # table(quant$peptide_id %in% psm$peptide_id)
  # quant %>% filter( ! peptide_id %in% psm$peptide_id)

  files_psm = fragpipe_find_psm_files(path)
  if(length(files_psm) == 0) {
    append_log(sprintf('cannot find any "psm.tsv" files in FragPipe results at "%s"', path), type = "error")
  }

  #### parse psm.tsv files
  psm = list()
  for(f in files_psm) {
    psm[[length(psm) + 1]] = fragpipe_parse_psm(f)
  }
  psm = bind_rows(psm)
  # enforce unique sample_id*peptide_id combinations (e.g. if multiple psm.tsv files contain data from rawfiles with the same name. Should not occur, but check anyway)
  psm = psm %>% distinct(sample_id, peptide_id, .keep_all = TRUE)
  # create "modified sequence" strings that are compatible with the MSstats.csv table so we can cross-reference
  psm = fragpipe_modseq_compose(psm)
  # quality control: check if any samples listed in the msstats table are missing from the psm.tsv tables
  tmp = setdiff(msstats_sampleid_unique, unique(psm$sample_id))
  if(length(tmp) > 0) {
    append_log(paste0("there is no data in any of the psm.tsv files for the following samples that are listed in MSstats.csv: ",
                      paste(tmp, collapse = ", "), "\nThis suggests your dataset is incomplete / some psm.tsv files are missing !"), type = "warning")
  }

  #### combine quant and psm data
  tib_result = left_join(quant, psm %>% select(sample_id, peptide_id, confidence, rt), by = c("sample_id", "peptide_id"))
  tib_result = fragpipe_map_combined_protein(tib_result %>% mutate(protein_id_long = protein_id), file_combined_protein)

  return(fragpipe_peptideresults_to_dataset(tib_result, acquisition_mode, confidence_threshold, collapse_peptide_by))
}



#' Helper function to find psm.tsv files within a FragPipe output directory
#'
#' @description
#' File "filelist_ionquant.txt" contains a list of all psm.tsv files that are part of this dataset.
#' In older FragPipe versions, e.g. 15.0, this file is not available.
#'
#' If filelist_ionquant.txt is available, extract listed psm.tsv files. If ALL of these are valid paths
#' for files within the provided FragPipe output directory (`path` parameter); assume this is the list
#' of PSM files we need. Print a warning if the table is in invalid format or if any of the listed
#' PSM files cannot be found within `path`.
#'
#' As a fallback, search all direct subdirectories (i.e. not recursive) of `path` for psm.tsv files.
#'
#' @param path FragPipe output directory
fragpipe_find_psm_files = function(path) {
  files_psm = NULL

  file_filelist_ionquant = path_exists(path, "filelist_ionquant.txt", try_compressed = TRUE, silent = TRUE)
  if(file_filelist_ionquant == "") { # fallback in case the file extension is updated in future versions
    file_filelist_ionquant = path_exists(path, "filelist_ionquant.tsv", try_compressed = TRUE, silent = TRUE)
  }

  if(file_filelist_ionquant != "") {
    # filelist_ionquant.txt suggests a plain text file, but it actually holds a tsv table with columns "flag" and "value"
    df = as_tibble(read_textfile_compressed(file_filelist_ionquant, as_table = TRUE))
    if(is.data.frame(df) && all(c("flag", "value") %in% colnames(df))) {
      if(any(df$flag == "--psm")) { # PSM files should have a special flag. If available, use it to subset
        df = df %>% filter(flag == "--psm")
      }
      tmp = grep("psm\\.tsv$", df$value, ignore.case = TRUE, value = TRUE)
      if(length(tmp) == 0) {
        append_log(sprintf('unexpected data found in "%s"; this should contain a list of psm.tsv files in the provided FragPipe output directory, however, no matching rows are found in this file/table. As a fallback, this function will search for psm.tsv within all subdirectories of the provided FragPipe output folder.', file_filelist_ionquant), type = "warning")
      }
      tmp_assume_relative = path_clean_slashes(paste0(path, "/", tmp))
      if(all(file.exists(tmp_assume_relative))) {
        files_psm = tmp_assume_relative
      } else {
        # never check the provided paths as-is, because those may be relative to the current working directory in the current R session
        # instead, as a fallback, assume the psm.tsv files are at most nested 1 dir-deep
        tmp_dir = dirname(tmp) # parent directory of provided psm.tsv file
        tmp_dir[tmp_dir == "" | tmp_dir == "/" | tmp_dir == "\\" | tmp_dir == "."] = "" # check for absence of parent dir, e.g. "psm.tsv" or "\\psm.tsv"
        tmp_assume_absolute = path_clean_slashes(paste0(path, "/", ifelse(tmp_dir == "", "", paste0(tmp_dir, "/")), basename(tmp)))
        if(all(file.exists(tmp_assume_absolute))) {
          files_psm = tmp_assume_absolute
        }
      }

      if(is.null(files_psm)) {
        append_log(sprintf('unexpected data found in "%s"; this should contain a list of psm.tsv files in the provided FragPipe output directory, however, some psm.tsv files listed in this table missing from the FragPipe output directory. As a fallback, this function will search for psm.tsv within all subdirectories of the provided FragPipe output folder.', file_filelist_ionquant), type = "warning")
      } else {
        append_log('extracted the list of psm.tsv files from "filelist_ionquant.txt"', type = "info")
      }
    } else {
      append_log(sprintf('unexpected data found in "%s"; this should be a table with columns "flag" and "value" and contain a list of psm.tsv files. As a fallback, this function will search for psm.tsv within all subdirectories of the provided FragPipe output folder.', file_filelist_ionquant), type = "warning")
    }
  }

  # Fallback 1: search all direct subdirectories (i.e. not recursive) of `path` for psm.tsv files.
  if(is.null(files_psm)) {
    subdirs = list.dirs(path, full.names = TRUE, recursive = FALSE)
    for(d in subdirs) {
      fpsm = path_exists(d, "psm.tsv", try_compressed = TRUE, silent = TRUE)
      if(fpsm != "") {
        files_psm = c(files_psm, fpsm)
      }
    }
    if(!is.null(files_psm)) {
      append_log(paste0(ifelse(file_filelist_ionquant == "", "", "fallback: "), "found psm.tsv files in subdirectories of the FragPipe output folder"), type = "info")
      # append_log(paste('Found psm.tsv files in subdirectories of the FragPipe output folder;', paste(basename(dirname(files_psm)), collapse=", ")), type = "info")
    }
  }

  # Fallback 2: check for a psm.tsv in the root directory (in typical fragpipe runs, that file should not exist, but check anyway)
  if(is.null(files_psm)) {
    tmp = path_exists(path, "psm.tsv", try_compressed = TRUE, silent = TRUE)
    if(tmp != "") {
      files_psm = tmp
      append_log("fallback: found psm.tsv file in the FragPipe output folder / root directory", type = "info")
    }
  }

  return(files_psm)
}



#' Alternative FragPipe workflow, operates directly on a single psm.tsv file
#'
#' @description
#' For typical FragPipe workflows, `import_dataset_fragpipe_ionquant()` should be used instead of this function !
#'
#' This function directly parses data from a psm.tsv file and doesn't use MBR data produced by FragPipe;
#' use this function for legacy FragPipe data  OR  FragPipe runs where no experiment info was provided at all in
#' the FragPipe workflow tab.
#'
#' if intensity_sum is set to FALSE, the most abundant PSM per peptide*sample is selected to represent abundance.
#'
#' @param filename the full file path of the input file
#' @param intensity_sum boolean value whether to take the sum intensity of all PSM per precursor (TRUE, default) or use the value from the PSM with highest abundance (FALSE)
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_fragpipe_psm_file = function(filename, intensity_sum = TRUE, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()
  append_log("importing FragPipe data from a single psm.tsv file using function import_dataset_fragpipe_psm_file(). If you only have 1 psm.tsv file from a FragPipe analysis without an MSstats.tsv file this is ok. But note that this function does not import IonQuant data, itensities are directly taken from FragPipe psm.tsv fiels!\nIn most cases you'll want to use MS-DAP function import_dataset_fragpipe_ionquant() instead, which imports all IonQuant data including MBR abundance values.", type = "warning")

  ### check for valid function parameters
  # must be a single string and point to an existing file (allow variants pointing to the compressed file)
  check_parameter_is_string(filename)
  if(dir.exists(filename)) {
    append_log('for this function, provide the full path to a single FragPipe psm.tsv (not a directory)', type = "error")
  }
  filename = path_exists(filename, NULL, try_compressed = TRUE, silent = FALSE) # throw error if not found

  check_parameter_is_string(acquisition_mode)
  if( ! acquisition_mode %in% c("dda","dia")) {
    append_log('acquisition_mode parameter must be either "dda" or "dia"', type = "error")
  }

  check_parameter_is_boolean(intensity_sum)
  check_parameter_is_numeric(confidence_threshold, minval_ = 0, maxval_ = 1)
  check_parameter_is_string(collapse_peptide_by)
  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }

  append_log("reading FragPipe PSM report...", type = "info")
  psm = fragpipe_parse_psm(filename, ifelse(intensity_sum, "min_confidence_sum_intensity", "min_confidence_max_intensity")) %>%
    select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, charge, confidence, rt, intensity) %>%
    filter(is.finite(intensity) & intensity > 0)

  return(fragpipe_peptideresults_to_dataset(psm, acquisition_mode, confidence_threshold, collapse_peptide_by))
}



#' Reconstruct modified peptide sequences from disjoint information provided in psm.tsv to facilitate cross-reference matching of between FragPipe output tables
#'
#' FragPipe documentation on modified sequences (that we here replace); https://github.com/Nesvilab/MSFragger/issues/138
#'
#' @description
#' Example data @ unchopped table in this function
#' HQGVMVGMGQK           5M(15.9949), 8M(15.9949) 5M(15.9949)  FALSE    FALSE    5     15.9949
#' HQGVMVGMGQK           5M(15.9949), 8M(15.9949) 8M(15.9949)  FALSE    FALSE    8     15.9949
#' ATEHPEPPK             N-term(42.0106)              N-term(42.0106)
#' ATTATMATSGSAR         6M(15.9949), N-term(42.0106) 6M(15.9949)
#' ATTATMATSGSAR         6M(15.9949), N-term(42.0106) N-term(42.0106)
#'
#' @param x psm table obtained via `fragpipe_parse_psm()`
fragpipe_modseq_compose = function(x) {
  # signature/ID to uniquely identify modified sequences
  x$temp_match_id = paste(x$sequence_plain, x$mods)

  x_result = x %>% filter(mods != "")
  # if there is nothing to do, quit  (sequence_modified should already be equal to sequence_plain @ upstream input validation)
  if(nrow(x_result) == 0) {
    return(x)
  }

  x_result = x_result %>%
    select(temp_match_id, sequence_plain, mods) %>%
    distinct(temp_match_id, .keep_all = TRUE) %>%
    mutate(
      uid = 1L:(dplyr::n()),
      sequence_modified_list = strsplit(sequence_plain, "", fixed = TRUE),
      mods_list = strsplit(mods, " *, *"), # comma delimited, but also deal with (accidental) surrounding whitespace
      first_aa = unlist(lapply(sequence_modified_list, head, n=1), recursive = FALSE, use.names = FALSE),
      last_aa = unlist(lapply(sequence_modified_list, tail, n=1), recursive = FALSE, use.names = FALSE)
    )

  x_unlisted = x_result %>%
    unchop(cols = mods_list) %>%
    mutate(
      is_nterm = grepl("N-term", mods_list, fixed = TRUE),
      is_cterm = grepl("C-term", mods_list, fixed = TRUE),
      index = sub("^(\\d+)\\D.*", "\\1", mods_list),
      new_aa_string = sub("(", "[", sub(")", "]", mods_list, fixed = TRUE), fixed = TRUE)
      # mass = sub(".*\\(([0-9.]+)\\).*", "\\1", mods_list)
    )

  x_unlisted$index[x_unlisted$is_nterm] = "1"
  x_unlisted$index[x_unlisted$is_cterm] = as.character(nchar(x_unlisted$index[x_unlisted$is_cterm]))
  x_unlisted$index = suppressWarnings(as.integer(x_unlisted$index))
  # z$mass = suppressWarnings(as.numeric(z$mass))
  rows_fail = is.na(x_unlisted$index) # | is.na(z$mass)
  if(any(rows_fail)) {
    append_log(sprintf("failed to extract sequence modification specs for %d unique sequences.\nExamples; %s", sum(rows_fail), paste(head(unique(x_unlisted$temp_match_id[rows_fail]), n=25), collapse = " , ")), type = "info")
  }

  x_unlisted$new_aa_string = sub(".*\\d([A-Z]\\[.*)", "\\1", x_unlisted$new_aa_string) # will fail / do nothing @ n-term and c-term mods
  x_unlisted$new_aa_string[x_unlisted$is_nterm] = paste0(sub("N-term", "n", x_unlisted$new_aa_string[x_unlisted$is_nterm]), x_unlisted$first_aa[x_unlisted$is_nterm])
  x_unlisted$new_aa_string[x_unlisted$is_cterm] = paste0(x_unlisted$last_aa[x_unlisted$is_nterm], sub("C-term", "c", x_unlisted$new_aa_string[x_unlisted$is_cterm]))
  # debug; z %>% filter(grepl("term", mods))

  for(i in 1L:nrow(x_unlisted)) {
    # directly update indices in the main, non-chopped, table @ index = x_unlisted$uid
    x_result$sequence_modified_list[[ x_unlisted$uid[i] ]][ x_unlisted$index[i] ] = x_unlisted$new_aa_string[i]
  }

  # from an aminoacid list to a string
  x_result$sequence_modified = unlist(lapply(x_result$sequence_modified_list, paste, collapse=""), recursive = FALSE, use.names = FALSE)

  # map modified sequences (main result) back to input table
  x$sequence_modified = x_result$sequence_modified[match(x$temp_match_id, x_result$temp_match_id)]
  rows = is.na(x$sequence_modified) # failed to match = not in modified sequence table = plain sequences
  x$sequence_modified[rows] = x$sequence_plain[rows]
  x$temp_match_id = NULL # remove temp column
  x$peptide_id = paste(x$sequence_modified, x$charge, sep = "_") # finally, update the composite peptide*charge ID
  return(x)
}



#' Parse a FragPipe MSstats.csv table
#'
#' @param filename full path to a FragPipe MSstats.csv file
fragpipe_parse_msstats = function(filename) {
  append_log(paste("parsing FragPipe file;", filename), type = "info")
  # read table and check whether all required columns are present
  x = read_table_by_header_spec(
    file = filename,
    attributes_required = list(
      protein_id = c("ProteinName", "Protein"),
      sequence_modified = c("PeptideSequence", "modified_peptide", "modified sequence"),
      charge = c("PrecursorCharge", "charge"),
      sample_id = "Run",
      intensity = "intensity"
    ),
    as_tibble_type = T
  ) %>%
    # enforce numeric column types and remove invalid entries (e.g. missing intensity values)
    mutate(
      charge = suppressWarnings(as.integer(charge)),
      intensity = suppressWarnings(as.numeric(intensity))
    ) %>%
    filter(!is.na(protein_id) & nchar(protein_id) > 3 &
             !is.na(sequence_modified) & nchar(sequence_modified) > 5 &
             is.finite(intensity) & intensity > 0)

  # remove mods, also remove lower-case symbol; "n[42.0106]STVHEILC[57.0215]K"
  useq = unique(x$sequence_modified)
  useq_plain = gsub("(\\[[^]]*\\])|(\\([^)]*\\))|[a-z]", "", useq) # cbind(useq, useq_plain)[useq != useq_plain,]
  x$sequence_plain = useq_plain[match(x$sequence_modified, useq)]
  x$peptide_id = paste(x$sequence_modified, x$charge, sep = "_")

  # enforce unique sample_id*peptide_id combinations (should not be needed, but check/enforce anyway)
  x %>% arrange(desc(intensity)) %>% distinct(sample_id, peptide_id, .keep_all = TRUE)
}



#' Parse a FragPipe quant.csv table
#'
#' @param filename full path to a FragPipe <rawfile>_quant.csv file
fragpipe_parse_quant = function(filename) {
  append_log(paste("parsing FragPipe file;", filename), type = "info")
  # read table and check whether all required columns are present
  x = read_table_by_header_spec(
    file = filename,
    attributes_required = list(
      sequence_plain = c("peptide", "peptide sequence", "sequence"),
      sequence_modified = c("modified_peptide", "modified sequence"),
      charge = "charge",
      rt = c("retention_time", "apex_retention_time"),
      intensity = "intensity"
    ),
    attributes_optional = list(
      mz = "calibrated_precursor_MZ"
    ),
    as_tibble_type = T
  )
  if(!"mz" %in% colnames(x)) {
    x$mz = NA_real_
  }

  x = x %>%
    # enforce numeric column types and remove invalid entries (e.g. missing intensity values)
    mutate(
      charge = suppressWarnings(as.integer(charge)),
      mz = suppressWarnings(as.numeric(mz)),
      rt = suppressWarnings(as.numeric(rt)),
      intensity = suppressWarnings(as.numeric(intensity))
    ) %>%
    filter(!is.na(sequence_plain) & nchar(sequence_plain) > 5 &
             is.finite(charge) & is.finite(rt) & is.finite(intensity)) %>%
    mutate(
      sequence_modified = ifelse(is.na(sequence_modified) | sequence_modified == "", sequence_plain, sequence_modified),
      peptide_id = paste(sequence_modified, charge, sep = "_")
    )

  # finally, sum intensities for the same peptide_id
  x %>%
    # highest intensity on top; for each peptide_id we'll return the m/z and RT matching this peak / precursor ion
    arrange(desc(intensity)) %>%
    group_by(peptide_id) %>%
    summarise(
      sequence_plain = sequence_plain[1],
      sequence_modified = sequence_modified[1],
      charge = charge[1],
      mz = mz[1],
      rt = rt[1],
      intensity = sum(intensity)
    ) %>%
    ungroup() %>%
    mutate(sample_id = sub("_quant.csv$", "", basename(filename), ignore.case = TRUE))
}



#' Parse a FragPipe psm.tsv table
#'
#' @param filename full path to a FragPipe psm.tsv file
#' @param precursor_unique_results how to deal with duplicate peptide_id*sample_id entries.
#' fast_min_confidence = just return min confidence (default: fast, but note that this ignores intensity values altogether)
#' min_confidence_max_intensity find the local minimum/maximum of respective values (pretty fast, bit slower than first option)
#' min_confidence_sum_intensity find the local minimum for confidence, summarize the intensities (slow)
fragpipe_parse_psm = function(filename, precursor_unique_results = "fast_min_confidence") {
  append_log(paste("parsing FragPipe file;", filename), type = "info")

  x = read_table_by_header_spec(
    file = filename,
    attributes_required = list(
      protein_id = "Protein",
      filename = "Spectrum File",
      sequence_plain = "Peptide",
      sequence_modified = "Modified Peptide",
      theoretical_mass = "Calculated Peptide Mass",
      mods = "Assigned Modifications",
      charge = "Charge",
      confidence = c("PeptideProphet.Probability", "Probability"),
      rt = "Retention",
      intensity = "Intensity"
    ),
    as_tibble_type = T
  ) %>%
    # enforce numeric column types and remove invalid entries (e.g. missing intensity values)
    mutate(
      charge = suppressWarnings(as.integer(charge)),
      intensity = suppressWarnings(as.numeric(intensity)),
      confidence = suppressWarnings(as.numeric(confidence))
    ) %>%
    filter(!is.na(protein_id) & nchar(protein_id) > 3 &
             !is.na(sequence_plain) & nchar(sequence_plain) > 5 &
             is.finite(charge) & is.finite(rt) & is.finite(confidence)) %>%
    mutate(
      sequence_modified = ifelse(is.na(sequence_modified) | sequence_modified == "", sequence_plain, sequence_modified),
      peptide_id = paste(sequence_modified, charge, sep = "_")
    )

  x$intensity[!is.finite(x$intensity)] = 0
  # In this R package we use confidence score = qvalue (or pvalue). So here take 1 - "PeptideProphet Probability"
  x$confidence[x$confidence < 0] = 0
  x$confidence[x$confidence > 1] = 1
  x$confidence = 1 - x$confidence

  # specifically for fragpipe; strip full path + "interact-" + ".pep.xml" extension
  # example entry in input table; C:\DATA\fragpipe_test\interact-a05191.pep.xml
  ufilename = unique(x$filename)
  ufilename_sampleid = sub("^interact[^a-z](.*)\\.pep\\.xml$", "\\1", basename(ufilename))
  x$sample_id = ufilename_sampleid[match(x$filename, ufilename)]


  x = x %>%
    mutate(tempid_sample_peptide = paste(sample_id, peptide_id)) %>%
    arrange(confidence, desc(intensity))

  # return 1 confidence value per sample*peptide_id
  # Note that this'll return 0 values for some intensities (e.g. non-zero value is available in another PSM with worse confidence)
  # when precursor_unique_results=="fast_min_confidence", we're done after this line
  result = x %>% distinct(tempid_sample_peptide, .keep_all = TRUE)

  if(precursor_unique_results == "min_confidence_max_intensity") {
    # result already has the best confidence per 'tempid_sample_peptide' identifier
    # now find the max intensity for this ID
    lookup = x %>% arrange(desc(intensity)) %>% distinct(tempid_sample_peptide, .keep_all = TRUE)
    result$intensity = lookup$intensity[match(result$tempid_sample_peptide, lookup$tempid_sample_peptide)]
    return(result %>% select(-tempid_sample_peptide))
  } else if(precursor_unique_results == "min_confidence_sum_intensity") {
    # result already has the best confidence per 'tempid_sample_peptide' identifier
    # now find the max intensity for this ID
    lookup = x %>%
      group_by(tempid_sample_peptide) %>%
      summarise(intensity = sum(intensity)) %>% # upstream we replaced NA with 0
      ungroup()
    result$intensity = lookup$intensity[match(result$tempid_sample_peptide, lookup$tempid_sample_peptide)]
  }

  return(result %>% select(-tempid_sample_peptide))
}



#' Read FragPipe combined_protein.tsv table and then update all protein_id in provided peptide table
#'
#' @param x peptide table composed at the end of `import_dataset_fragpipe_ionquant()` or `import_dataset_fragpipe_quantfiles()`
#' @param file_combined_protein see `import_dataset_fragpipe_ionquant()`
fragpipe_map_combined_protein = function(x, file_combined_protein) {

  ### updated protein_id with ambiguous included
  protid_to_fullprotid = parse_fragpipe_proteins_from_combined_proteins(file_combined_protein)
  # ensure the protein table covers all proteins in our precursor-level result table
  pid_mismatch = setdiff(unique(x$protein_id_long), protid_to_fullprotid$protein_id_long)
  if(length(pid_mismatch) > 0) {
    append_log(sprintf('%d unique protein ID from psm.tsv and ion.tsv files could not be matched against combined_protein.tsv , so "Indistinguishable Proteins" are not available for these protein groups (does not affect MS-DAP data processing, only affects assigned protein and gene ID for respective protein groups which now only describe the leading protein)', length(pid_mismatch)), type = "warning")
  }

  # update protein_id in peptide table
  x %>%
    left_join(protid_to_fullprotid, by = "protein_id_long") %>%
    # if we found the protein_id in the master table, update protein-group-ID
    mutate(protein_id = ifelse(is.na(protein_id_long_plus_ambiguous), protein_id_long, protein_id_long_plus_ambiguous)) %>%
    select(-protein_id_long, -protein_id_long_plus_ambiguous)
}



#' From FragPipe peptide result table to an MS-DAP dataset
#'
#' @description
#' After composing a protein_id,peptide_id,sample_id,intensity,etc... table upstream, fix protein identifiers,
#' transform RT from seconds to minutes, transform intensity to log2, collapses peptides by plainseq/modseq and finally construct a dataset.
#'
#' @param tib_result peptide table composed at the end of `import_dataset_fragpipe_ionquant()` or `import_dataset_fragpipe_quantfiles()`
#' @param acquisition_mode see `import_dataset_fragpipe_ionquant()`
#' @param confidence_threshold see `import_dataset_fragpipe_ionquant()`
#' @param collapse_peptide_by see `import_dataset_fragpipe_ionquant()`
fragpipe_peptideresults_to_dataset = function(tib_result, acquisition_mode, confidence_threshold, collapse_peptide_by) {

  # unpack the protein identifiers > derive short IDs > flag decoys > summarize at proteingroup level
  tib_uprot = proteingroup_to_idshort_lookuptable(tib_result$protein_id)
  # overwrite result table protein identifiers
  i = match(tib_result$protein_id, tib_uprot$protein_id)
  tib_result$protein_id = tib_uprot$idshort[i]
  tib_result$isdecoy = tib_uprot$isdecoy[i]


  ### format peptide table columns, e.g. ensure intenties are log-transformed
  tib_result = tib_result %>%
    # remove invalid protein_id (empty/NA, as flagged by idshort lookup table  AND  remove decoys)
    filter(!is.na(protein_id) & isdecoy == FALSE) %>%
    mutate(
      detect = is.finite(confidence) & confidence <= confidence_threshold, # define 'detect'
      rt = rt / 60, # convert RT to minutes
      intensity = log2(intensity) # log2 transform intensities
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

  x %>%
    select(protein_id_long = Protein, protein_id_long_ambiguous = `Indistinguishable Proteins`) %>%
    distinct(protein_id_long, .keep_all = TRUE) %>% # shouldn't be needed, but enforce unique protein IDs in the table we just read from file
    mutate(
      protein_id_long_plus_ambiguous = ifelse(is.na(protein_id_long_ambiguous) | nchar(protein_id_long_ambiguous) < 5,
                                              protein_id_long,
                                              paste(protein_id_long, protein_id_long_ambiguous, sep = ";")),
      # standardize separation characters, but only for the "replacement" identifiers (and not the IDs that are to be matched to other FragPipe tables)
      protein_id_long_plus_ambiguous = gsub("[ ,;]+", ";", protein_id_long_plus_ambiguous)
    ) %>%
    # column names; make explicit that the protein identifiers returned by this function are in long-format, i.e. not the default in MS-DAP
    select(protein_id_long, protein_id_long_plus_ambiguous)
}
