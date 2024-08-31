
#' Import sample metadata from an Excel table
#'
#' @param dataset your dataset
#' @param filename full file path of the input file (eg; C:/temp/template_experiment1.xlsx)
#'
#' @export
import_sample_metadata = function(dataset, filename) {
  # TODO: input validation
  stopifnot(is.list(dataset) && "peptides" %in% names(dataset))

  # reset cache
  dataset = invalidate_cache(dataset)

  dataset$samples = sample_metadata_from_file(sample_id = unique(dataset$peptides$sample_id), filename = filename)
  # add peptide and protein counts to samples tibble, so downstream code includes this in output tables and report figures relating to sample metadata
  dataset$samples = peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))
  return(dataset)
}



#' For advanced users: import sample metadata directly from the sample names using regular expressions
#'
#' @param dataset your dataset
#' @param sample_property_regex a named list of regular expressions to extract sample properties. name=sample property, value = array of length two with the regex and the replace string (first 2 arguments of sub())
#' @param sample_exclude_regex optionally, regex to select samples that should be flagged as 'exclude'
#' @param group_order optionally, you can provide a preferred ordering of the sample groups as a character array
#' @param group_regex_array optionally, a regex that specifically extracts sample groups and leaves non-matches in an 'unassigned' group
#'
#' @examples
#' # example of a regex applied to sample names to extract the sample shortname and group
#' \dontrun{sample_property_regex = list(
#'   shortname = c(".*experimentname_([A-Z]+.\\d).*", "\\1"),
#'   group = c(".*experimentname_([A-Z]+)_.*", "\\1")
#'   )
#' }
#' @export
sample_metadata_custom = function(dataset, sample_property_regex = list(), sample_exclude_regex = "", group_order = NA, group_regex_array = NA) {
  # TODO: input validation
  stopifnot(is.list(dataset) && "peptides" %in% names(dataset))

  dataset$samples = sample_metadata_from_filenames(sample_id = unique(dataset$peptides$sample_id), sample_property_regex=sample_property_regex, sample_exclude_regex=sample_exclude_regex, group_order=group_order, group_regex_array=group_regex_array)
  # add peptide and protein counts to samples tibble, so downstream code includes this in output tables and report figures relating to sample metadata
  dataset$samples = peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))
  return(dataset)
}



#' Write an Excel sheet to file that contains all samples in the dataset
#'
#' @param dataset your dataset
#' @param filename full file path of the desired output file (eg; C:/temp/template_experiment1.xlsx)
#' @param overwrite logical indicating whether to overwrite existing file (if any)
#'
#' @importFrom gtools mixedsort
#' @importFrom openxlsx write.xlsx
#' @export
write_template_for_sample_metadata = function(dataset, filename, overwrite = FALSE) {
  if(!grepl("\\.xlsx$", filename, ignore.case = T)) {
    filename = paste0(filename, ".xlsx")
    append_log(paste("sample metadata filename should end with .xlsx, appending extension. Updated filename:", filename), type = "warning")
  }

  if(nchar(filename) > 260) {
    append_log(paste("failed to create output file, total path length longer than 260 characters. Choose an output_dir with shorter path length, eg; C:/data/myproject (Windows) or /home/user/myproject (unix). Attempted output file (longer than 260 characters, so invalid);", filename), type="error")
  }

  if(file.exists(filename)) {
    if(overwrite) {
      remove_file_if_exists(filename)
      append_log(paste("overwriting file:", filename), type = "warning")
    } else {
      append_log(sprintf("file already exists: %s. Either remove file, or call this function with overwrite=TRUE parameter", filename), type = "error")
    }
  } else {
    parent_dir = dirname(filename)
    if(!dir.exists(parent_dir)) {
      dir.create(parent_dir, recursive = T)
    }
  }

  # get the unique set of sample IDs and apply natural sort
  s = gtools::mixedsort(unique(dataset$peptides$sample_id))
  # remove leading/trailing special characters  +  strip substrings at start or end that are found in every sample_id
  # tolower() makes this code much more robust, in practise we find lots of "dirty" input where filenames inconsistently mix casing (eg; WT1 wt2 WT3 ko1)
  s_clean = gsub("(^[^0-9a-z]+)|([^0-9a-z]+$)", "", tolower(s)) # example; s = c("WT1_", "WT_2", ".WT3", "WT3-")
  s_clean = strip_common_substring(s_clean) # example: s_clean = c("20191231_WT1_bob", "20191231_WT2_bob", "20191231_KO1_bob", "20191231_KO2_bob")
  s_clean = gsub("(^[^0-9a-z]+)|([^0-9a-z]+$)", "", s_clean) # example; "8" "_1" "_2" "_3"
  # shortname; replace special characters with underscore to encourage clean names. most common sep char -->> rename to underscore, all other to "-"
  s_shortname = s_clean # lowercase converted upstream, so regex' can be case sensitive
  sepchar_most_common = count_sep_char(s_shortname)
  if(nrow(sepchar_most_common) > 0) {
    s_shortname = gsub(sepchar_most_common$char_regex[1], "$", s_shortname) # most common sep char -->> $
    s_shortname = gsub("[^0-9a-z$]", "-", s_shortname) # every special char except $ -->> -
    s_shortname = gsub("$", "_", s_shortname, fixed = T) # $ -->> _
  }

  # sample metadata table
  df = data.frame(sample_id=s, shortname=s_shortname, exclude="FALSE", group=rep("", length(s)), stringsAsFactors = F)
  # infer / guess experiment metadata encoded in sample names. Expect "dirty" input that is not consistently formatted
  # start with s_clean so different types of special characters are maintained, this function will use it to infer if user gave different meaning to distinct separation chars (eg; maybe whitespace is separation character and underscore is not. "condition_a rep2" "condition2 rep2")
  s_metadata = metadata_matrix_from_filenames(s_clean)
  if(is.matrix(s_metadata) && ncol(s_metadata) > 1) {
    # if all is well, s_metadata contains nicely formatted/split metadata but some special chars may remain -->> regex replace (allow - though, maybe these are formatted dates)
    s_metadata = apply(s_metadata, 2, function(x) gsub("[^0-9a-z_-]+", "_", x, ignore.case = T))
    df = cbind(df, s_metadata) # converts to matrix, but that's fine
  }
  # openxlsx::write.xlsx(df, file = filename)  # simple case, just write data table (forego making a formatted Excel doc that includes help/documentation)

  wb = openxlsx::createWorkbook()
  header_style = openxlsx::createStyle(textDecoration = "Bold")
  openxlsx::addWorksheet(wb, "samples")
  openxlsx::addWorksheet(wb, "instructions")
  openxlsx::writeData(wb, "samples", df, keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "instructions", data.frame(
    Instructions = c("This table should describe metadata for all samples in your dataset.",
                     "",
                     "sample_id (required): lists the input sample names, auto-generated from your dataset, should not be edited.",
                     "shortname (required): a recognizable name/label for each sample that will be used in QC plots. MS-DAP tried to infer this, please curate. Allowed characters: alphanumeric, underscore, space, plus and minus. We encourage using simple/clean names with just alphanumeric and underscores.",
                     "fraction (optional):  column for datasets with fractionated data. If you add this column, data from samples with the same shortname are merged. Do not add column with this name if your data is NOT fractionated ! empty values not allowed.",
                     "group (required): describes the sample groups, these labels should group samples at similar experimental conditions. Allowed characters: alphanumeric, underscore or minus.",
                     "exclude (optional): boolean flag that indicates which samples you consider to be strong outliers that should be excluded from statistical analyses (but will be included in QC figures, highlighted as 'exclude' status)",
                     "",
                     "QC figures will be automatically generated by MS-DAP for all sample metadata you provide !",
                     "So we highly recommend describing additional sample metadata by simply adding more columns beyond those listed above (only for parameters/columns that describe more than 1 unique value of course).",
                     "Examples of metadata to consider adding; experiment date or batch, order of sample handling, order measurement, gel number, gel lanes, cell counts for FACS experiment, etc.",
                     "Allowed characters are: alphanumeric, underscore, space, dot, minus, (forward)slash, backslash. We recommend sticking with alphanumeric, underscore and minus characters.",
                     "",
                     "MS-DAP tries to guess/infer metadata from sample names to convenience data entry. If successful, these were added in columns that start with an 'x'.",
                     "Please double-check and rename respective columns to meaningful names/labels.",
                     "In most cases, one represents your sample groups (move data to 'group' column and remove respective 'x' column), one the replicate IDs (we typically remove these, redundant info with 'shortname' column), and others remaining experiment conditions (rename columns and curate !)",
                     "",
                     "reserved 'keywords'; column names should not start with 'condition' or 'contrast', nor have the name 'sample_index'")), keepNA = FALSE, headerStyle = header_style)
  openxlsx::saveWorkbook(wb, file = filename)

  append_log(sprintf('Sample metadata file "%s" has been created. Please open it with Excel/LibreOffice and enter experiment metadata (this file contains documentation at the "instructions" tab)', filename), type = "info")
}



#' extract metadata from sample names using regular expressions
#'
#' note; user-facing function is sample_metadata_custom()
#'
#' @param sample_id character array representing sample IDs
#' @param sample_property_regex a named list of regular expressions to extract sample properties. name=sample property, value = array of length two with the regex and the replace string (first 2 arguments of sub())
#' @param sample_exclude_regex optionally, regex to select samples that should be flagged as 'exclude'
#' @param group_order optionally, you can provide a preferred ordering of the sample groups as a character array
#' @param group_regex_array optionally, a regex that specifically extracts sample groups and leaves non-matches in an 'unassigned' group
sample_metadata_from_filenames = function(sample_id, sample_property_regex = list(), sample_exclude_regex = "", group_order = NA, group_regex_array = NA) {
  # input validation; either sample_property_regex or group_regex_array
  if(length(sample_property_regex) == 0 || !is.list(sample_property_regex)) {
    if (length(group_regex_array) == 0 || any(is.na(group_regex_array))) {
      append_log("provide either 'group_regex_array' (regex array, no NA's) or 'sample_property_regex' (regex list, no NA's) argument", type = "error")
    }
  } else {
    if (!all(c("shortname", "group") %in% names(sample_property_regex))) {
      append_log("'sample_property_regex' must be a list, required elements are 'shortname' and 'group'", type = "error")
    }
  }

  if ("sample_id" %in% names(sample_property_regex)) {
    append_log("you cannot re-define the 'sample_id' @ 'sample_property_regex'", type = "error")
  }

  sample_id = unique(sample_id)

  # unique sample names from peptide tibble
  df = data.frame(sample_id=sample_id, shortname=sample_id, stringsAsFactors = F)
  # add all attributes from regex
  for (n in names(sample_property_regex)) {
    n_regex = sample_property_regex[[n]]
    n_replacement = ""
    if (length(n_regex) > 1) {
      n_replacement = n_regex[2]
    }
    df[, n] = gsub(n_regex[1], n_replacement, df$sample_id, ignore.case = T)
  }

  if (length(group_regex_array) > 0 && !any(is.na(group_regex_array))) {
    df$group = "" # default; no groups
    # iterate regex array to find respective samples
    for (i in seq_along(group_regex_array)) {
      grp = names(group_regex_array)[i]
      rows = grepl(group_regex_array[i], df$sample_id, ignore.case = T)
      if (any(rows)) df$group[rows] = grp
    }
    if (any(df$group == "")) {
      df$exclude = df$group == ""
      df$group[df$exclude] = "_unassigned_"
      # print(df)
      # append_log("the group regex array is invalid, some samples have no associated group", type = "error")
    }
  }

  # QC & sort
  return(sample_metadata_sort_and_filter(df, sample_exclude_regex = sample_exclude_regex, group_order = group_order))
}



#' read sample metadata table from file and match it against provided sample IDs
#'
#' note; user-facing function is import_sample_metadata()
#'
#' @param sample_id character array representing sample IDs
#' @param filename full file path of the input file (eg; C:/temp/template_experiment1.xlsx)
#' @param group_order optionally, you can provide a preferred ordering of the sample groups as a character array
#'
#' @importFrom data.table fread
#' @importFrom openxlsx read.xlsx
sample_metadata_from_file = function(sample_id, filename, group_order = NA) {
  check_parameter_is_string(filename)
  if (!file.exists(filename)) {
    append_log(paste("File does not exist:", filename), type = "error")
  }
  is_type_tsv = grepl("\\.(csv|tsv)$", filename, ignore.case = T)
  is_type_xlsx = grepl("\\.(xlsx)$", filename, ignore.case = T)
  if (!is_type_tsv && !is_type_xlsx) {
    append_log(paste("Sample metadata should either be a csv/tsv file, or an Excel table stored as xlsx:", filename), type = "error")
  }

  ### read from disk
  if(is_type_xlsx) {
    df = openxlsx::read.xlsx(filename, sheet=1, sep.names="_", detectDates=T)
  } else {
    df = as.data.frame(data.table::fread(filename))
  }
  if (!is.data.frame(df) || nrow(df) == 0) {
    append_log("provided sample metadata file was empty", type = "error")
  }
  if (! "sample_id" %in% colnames(df)) {
    append_log("sample_id is a required attribute for sample metadata", type = "error")
  }

  ### before thorough downstream checking, basic QC check for code in this function; character and not NA
  tmp = as.character(df$sample_id)
  tmp[is.na(df$sample_id)] = ""
  df$sample_id = tmp
  # remove rows where sample_id column is empty
  df = df[df$sample_id != "NA" & df$sample_id != "", ]

  # check dupes for non-empty  (e.g. deal with trailing empty rows etc. in downstream code/functions)
  df_sid_dupe = df$sample_id != "" & duplicated(df$sample_id)
  if(any(df_sid_dupe)) {
    append_log(paste0("sample metadata table must not contain any duplicates in sample_id column. Duplicated entries;\n", paste(unique(df$sample_id[df_sid_dupe]), collapse="\n")), type = "error")
  }


  ### increase cross-compatibility so users don't have to edit samples.xlsx
  # maybe the table we read from file contains file extensions that we stripped out in (newer versions of) MS-DAP
  # so e.g. "sample1.wiff" in `df$sample_id` is present in `sample_id` parameter as "sample1"
  # solution: update `df` for elements that do not match currently, and have exactly 1 match after regex
  # strip path and whitelisted extensions from filename
  df_sid_strip = gsub(regex_rawfile_strip_extension(), "", df$sample_id, ignore.case = T)
  sample_id__matched = sample_id %in% df$sample_id
  for(i in which(df$sample_id != "")) {
    if(!any(df$sample_id[i] == sample_id)) { # df element i is unmatched
      index_in_param = which(sample_id == df_sid_strip[i]) # try i's stripped filename
      if(length(index_in_param) == 1 && !sample_id__matched[index_in_param]) { # exactly 1 match, and that element was not matched already
        df$sample_id[i] = df_sid_strip[i] # update user's metadata table from file, replacing the sample_id with stripped variant
        sample_id__matched[index_in_param] = TRUE # flag as matched
      }
    }
  }


  ### input table should only contain samples matching this dataset
  samples_missing = setdiff(sample_id, df[,"sample_id"])
  samples_extra = setdiff(df[,"sample_id"], sample_id)
  if (length(samples_missing) > 0) {
    append_log(paste0("sample metadata table must contain all samples from this dataset. Missing samples:\n", paste(samples_missing, collapse="\n")), type = "error")
  }
  if (length(samples_extra) > 0) {
    append_log(paste0("sample metadata table must contain only samples from this dataset. Additional samples that should be removed from metadata table:\n", paste(samples_extra, collapse="\n")), type = "error")
  }

  if(length(group_order) == 0 || any(!is.na(group_order))) {
    group_order = NA
  } else {
    group_order = unique(group_order)
  }

  # QC & sort
  return(sample_metadata_sort_and_filter(df, group_order = group_order))
}



#' sample metadata formatting and validation
#'
#' @param df data.frame with sample metadata
#' @param sample_exclude_regex optionally, regex to select samples that should be flagged as 'exclude'
#' @param group_order optionally, you can provide a preferred ordering of the sample groups as a character array
#'
#' @importFrom gtools mixedorder
#' @importFrom stringr str_squish
#' @importFrom tibble rowid_to_column
sample_metadata_sort_and_filter = function(df, sample_exclude_regex = "", group_order = NA) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    append_log("sample metadata must be a non-empty dataframe", type = "error")
  }

  ### first, format input table
  # convert all properties (/ column names) to lowercase  +  replace all non-alphanumerics with underscore
  colnames(df) = gsub("[^0-9a-z]+", "_", tolower(colnames(df)))
  colnames(df) = gsub("(^_)|(_$)", "", colnames(df)) # strip leading/trailing underscores

  # dupes in column names
  if (anyDuplicated(colnames(df)) != 0) {
    append_log(paste0("duplicated column names are not allowed. Duplicated column names (after conversion to lowercase and substituting special chars with underscores);\n", paste(colnames(df)[duplicated(colnames(df))], collapse = ", ")), type = "error")
  }

  # input validation on column names: contrasts must be constructed in-code (for now)
  cols_invalid = grepl("^condition", colnames(df), ignore.case = T)
  if (any(cols_invalid)) {
    append_log("metadata column names cannot start with \"condition\" (reserved keyword)", type = "error")
  }

  # input validation on column names: contrasts must be constructed in-code (for now)
  cols_remove = grepl("^contrast", colnames(df), ignore.case = T)
  if (any(cols_remove)) {
    df = df[ , !cols_remove, drop=F]
    append_log("importing contrasts from file is not supported yet, all columns starting with \"contrast\" have been removed", type = "warning")
  }

  # remove columns generated by ms-dap (eg; when re-importing a samples.xlsx previously generated by this pipeline)
  cols_remove = c(intersect("sample_index", colnames(df)), grep("^(detected|all)_(peptides|proteins)$", colnames(df), value=T, ignore.case = T))
  if (length(cols_remove) > 0) {
    df = df[ , ! colnames(df) %in% cols_remove, drop=F]
    append_log(paste("removed input columns (reserved keywords);", paste(cols_remove, collapse=", ")), type = "warning")
  }

  # trim whitespace, to prevent mistakes in input data (eg; group name has a trailing space in one row -->> should not be considered a distinct group)
  df = df %>% mutate(dplyr::across(where(is_character), stringr::str_squish)) # apply function to character columns; https://dplyr.tidyverse.org/reference/across.html

  # replace empty strings with NA. Importantly, do this after above whitespace trim so any cell with only spaces is recognised as empty
  df[is.na(df) | df == ""] = NA

  # silently drop completely empty columns (eg; Excel may generate 'invisible' extra columns padded at the end)
  # and also those with just 1 value (useless for sample metadata plots anyway)
  col_all_empty = colSums(!is.na(df)) <= as.integer(nrow(df) > 1) # if only 1 row in df, only drop empty columns. otherwise, drop columns that have 1 or 0 values
  if(any(col_all_empty)) {
    # append_log(sprintf("removing empty columns from sample metadata table;", paste(colnames(df)[col_all_empty], collapse=", ")), type = "warning")
    df = df[, !col_all_empty, drop=F]
  }


  ### input validation on contents of the table
  # required properties/columns
  if (!all(c("sample_id", "shortname", "group") %in% colnames(df))) {
    append_log("sample_id, shortname and group are required attributes for sample metadata", type = "error")
  }

  if (any(is.na(df[, c("sample_id", "shortname", "group")]))) {
    append_log("no empty values allowed in columns; sample_id, shortname and group", type = "error")
  }

  # dupes in sample_id
  if (anyDuplicated(df$sample_id) != 0) {
    print(df[duplicated(df$sample_id),])
    append_log("duplicate entries in the 'sample_id' column are not allowed", type = "error")
  }

  # input validation on sample_id: must not be integers. This protects against possibly common R implementation bugs; indexing by column names is prone to bugs when the column names are string representations of integers
  sample_id_invalid = grep("^\\d+$", df$sample_id, value=T)
  if (length(sample_id_invalid) > 0) {
    append_log(paste("We advice against numeric sample_id values. Preferably, use an alphanumeric string to prevent confusion between sample index and input filenames. Affected sample_id values;", paste(unique(sample_id_invalid), collapse = ", ")), type = "warning")
  }

  # analogous for shortname
  shortname_invalid = grep("^\\d+$", df$shortname, value=T)
  if (length(shortname_invalid) > 0) {
    append_log(paste("We advice against numeric shortname values. Preferably, use an alphanumeric string to prevent confusion between sample index and sample names in visualizations. Affected shortname values;", paste(unique(shortname_invalid), collapse = ", ")), type = "warning")
  }

  if ("fractions" %in% colnames(df) && !"fraction" %in% colnames(df)) {
    append_log("a common typo in sample metadata column names is using 'fractions' instead of 'fraction', fixing now by renaming this column...", type = "warning")
    colnames(df)[colnames(df)=="fractions"] = "fraction"
  }

  # check for empty fractions
  if ("fraction" %in% colnames(df) && any(is.na(df$fraction))) {
    append_log("empty values in the fraction column are not allowed", type = "error")
  }

  # check for dupes in fraction * shortname
  if ("fraction" %in% colnames(df) && anyDuplicated(df[, c("shortname", "fraction")]) != 0) {
    print(df[duplicated(df[, c("shortname", "fraction")]),])
    append_log("duplicated fraction * sample-name ('shortname' column) combinations are not allowed", type = "error")
  }

  if (!"fraction" %in% colnames(df) && anyDuplicated(df$shortname) != 0) {
    print(df[duplicated(df$shortname),])
    append_log("duplicated sample names in the 'shortname' column are not allowed (unless your samples were fractionated, in which case you should provide a column named 'fraction')", type = "error")
  }

  # replace repeated whitespace with single
  df = df %>% mutate_if(is.character, function(x) gsub("\\s+", " ", x))

  # enforce valid characters in metadata
  df_check = df[, setdiff(colnames(df), "sample_id")]
  props_invalid = colnames(df_check)[apply(df_check, 2, function(x) any(!is.na(x) & is.character(x) & grepl("[^0-9a-z_+., /;:#!?*-]", x, ignore.case = T) ) )] # apply regex to each column, check if any row matches; allow alphanumeric, underscore, space, minus, slash, dot
  if (length(props_invalid) > 0) {
    print(df[ , props_invalid, drop = F])
    append_log(paste0("invalid characters in sample metadata columns; ", paste(props_invalid, collapse = ", "), ". valid characters allowed are: alphanumeric, underscore, space (but no double spaces), dot, minus, (forward)slash, backslash, semicolon, hashtag, star"), type = "error")
  }
  rows_invalid = grepl("[^0-9a-z_+. -]", df$group, ignore.case = T) | grepl(" {2,}", df$group)
  if (any(rows_invalid)) {
    append_log(paste0("invalid characters in sample metadata, column group, values; ", paste0(unique(df$group[rows_invalid]), collapse=", "), " -->> group column may only contain characters that are; alphanumeric, underscore, space (but no double spaces), dot, plus and minus"), type = "error")
  }
  rows_invalid = grepl("[^0-9a-z_+. -]", df$shortname, ignore.case = T) | grepl(" {2,}", df$shortname)
  if (any(rows_invalid)) {
    append_log(paste0("invalid characters in sample metadata, column shortname, values; ", paste0(unique(df$shortname[rows_invalid]), collapse=", "), " -->> Shortnames may only contain characters that are; alphanumeric, underscore, space (but no double spaces), dot, plus and minus"), type = "error")
  }


  ### enforce data types: shortname and group must be character type
  # Important, otherwise all downstream code will need to be checked for silent !! type conversion errors
  # example;
  # some code matches columns from a peptide abundance matrix (matrix type, peptide_id*sample_id), character type, against dataset$samples$shortname to find respective metadata
  # -->> if shortname is numeric type, this yields no match but also no error (dynamic typing, R has no strong types, etc.)
  df$sample_id = as.character(df$sample_id)
  df$shortname = as.character(df$shortname)
  df$group = as.character(df$group)


  ### sort data
  # sort samples by shortname, then by group (use gtools::mixedorder to nicely handle not-super-cleanly formatted samples. eg; 10_WT, 3_WT, 05_WT)
  df = df[ gtools::mixedorder(df$shortname), , drop = F]
  df = df[ gtools::mixedorder(df$group), , drop = F]

  # next, sort by group if a preference for (a subset of) groups is given, sort by matching against this vector. eg; group_order=c("WT","KO")
  if (length(group_order) > 0 && all(!is.na(group_order))) {
    df = df[ order(match(tolower(df$group), tolower(group_order)), na.last = T), ]
  }


  ### exclude samples
  # flag unwanted files
  if(!"exclude" %in% colnames(df)) {
    # not in df yet, init as default false
    df$exclude = FALSE
  } else {
    # tolower is also a string conversion (so works with both character and numeric input). %in% deals with NA's
    df$exclude = tolower(df$exclude) %in% c("true", "1")
  }

  if (length(sample_exclude_regex) == 1 && !is.na(sample_exclude_regex) && sample_exclude_regex != "") {
    df$exclude = grepl(sample_exclude_regex, df$sample_id, ignore.case = T)
    if (!any(df$exclude)) {
      append_log("your regex for defining 'exclude' samples did not match anything", type = "warning")
    }
  }

  # there must be at least 1 sample group with 1 non-exclude sample
  if(all(df$exclude == TRUE)) {
    append_log("there must be at least 1 sample that is not flagged as 'exclude'", type = "error")
  }

  # finally, add index as the first column
  df = tibble::rowid_to_column(df, var = "sample_index")

  return(as_tibble(df))
}



#' test if the dataset has contrast definitions that are incompatible with MS-DAP version 1.2
#'
#' @param dataset your dataset
#' @export
has_legacy_contrast_definitions = function(dataset) {
  is.list(dataset) && (
    ("samples" %in% names(dataset) && is.data.frame(dataset$samples) && any(grepl("^contrast:", colnames(dataset$samples), ignore.case = TRUE))) ||
      # regex: has intensity_contrast column that does not contain a hashtag
      # (all new-style contrast names do because the dataset$samples column with "condition" is encoded in contrast labels since MS-DAP release 1.2)
      ("peptides" %in% names(dataset) && is.data.frame(dataset$peptides) && any(grepl("^intensity_contrast[^#]+$", colnames(dataset$samples), ignore.case = TRUE)))
  )
}



#' throw an error if has_legacy_contrast_definitions() is TRUE
#'
#' @param dataset your dataset
error_legacy_contrast_definitions = function(dataset) {
  if(has_legacy_contrast_definitions(dataset)) {
    append_log("this dataset was constructed with an older, incompatible MS-DAP version. To fix this, remove and then redefine all contrasts with the following 2 steps:\ndataset = remove_contrasts(dataset)\ndataset = setup_contrasts(dataset, ... same definitions you used before ...)\nDone! you may now use the updated dataset in MS-DAP as per usual", type = "error")
  }
}



#' Remove all contrast definitions and respective filtering and DEA results
#'
#' @description
#' also removes contrast-related information from the dataset that was used in older MS-DAP version (prior to version 1.2)
#'
#' @param dataset your dataset
#' @seealso [print_contrasts()], [add_contrast()]
#' @export
remove_contrasts = function(dataset) {
  if(!is.list(dataset)) {
    append_log("dataset parameter must be a list", type = "error")
  }
  # remove legacy random variables specification if it exists
  dataset$dea_random_variables = NULL
  # remove legacy contrast definitions in dataset$samples
  if("samples" %in% names(dataset) && is.data.frame(dataset$samples)) {
    dataset$samples = dataset$samples %>% select(-tidyselect::starts_with("contrast:"))
  }
  # remove all contrast-filters from peptide table
  if("peptides" %in% names(dataset) && is.data.frame(dataset$peptides)) {
    dataset$peptides = dataset$peptides %>% select(-tidyselect::starts_with("intensity_contrast"))
  }
  # remove all dea and dd results
  dataset$de_proteins = NULL
  dataset$dd_proteins = NULL
  # init empty list
  dataset$contrasts = list()
  return(dataset)
}



#' Print an overview of all defined contrasts and their respective sample metadata tables to the console
#'
#' @param dataset your dataset
#' @seealso [add_contrast()], [remove_contrasts()]
#' @export
print_contrasts = function(dataset) {
  if(!is.list(dataset)) {
    append_log("dataset parameter must be a list", type = "error")
  }

  if(!is.list(dataset$contrasts) || length(dataset$contrasts) == 0) {
    append_log("no contrasts have been defined", type = "info")
  }

  for(i in seq_along(dataset$contrasts)) {
    contr = dataset$contrasts[[i]]
    cat(ifelse(i > 1, "\n", ""), contr$label, "\n")
    print(inner_join(
      dataset$samples %>% select(sample_id, shortname), # optionally, include 'fraction': tidyselect::any_of("fraction")
      contr$sample_table,
      by = "sample_id"
    ))
  }
}



#' Create a contrast for differential expression analysis
#'
#' @description
#' Note that a MS-DAP contrast for "A vs B" will return foldchanges for B/A in downstream output
#' tables and data visualizations. For example, for the contrast "control vs disease" a positive
#' log2 foldchange implies protein abundances are higher in the "disease" sample group.
#'
#' Throughout this function, all samples where the column "exclude" is set to TRUE are disregarded.
#'
#' @examples \dontrun{
#'   # first, remove all existing contrasts (and their respecive DEA results)
#'   dataset = remove_contrasts(dataset)
#'
#'   # Assume that column "group" in your sample metadata table specifies sample groups
#'   # The following example will create the contrast "WT vs KO"
#'   dataset = add_contrast(dataset, "group", "WT", "KO")
#'
#'   # If the sample metadata table contains a column "batch" that should be used as
#'   # a regression variable (works for ebayes/deqms/msqrob), we can add it as follows:
#'   dataset = add_contrast(dataset, "group", "WT", "KO", colname_additional_variables = "batch")
#'
#'   # Elaborate example; create a contrast while matching multiple values per group.
#'   # Testing all motor- and visual-cortex samples (described in the "brain_region"
#'   # column) against prefrontal cortex samples, with additional regression variables
#'   # batch and age.
#'   # In this example, the sample metadata table must contain columns
#'   # "brain_region", "batch", "age"
#'   dataset = add_contrast(
#'     dataset,
#'     # this parameter describes 1 column name in `dataset$samples`
#'     colname_condition_variable = "brain_region",
#'     # values in "brain_region" that are the first group in A/B testing
#'     values_condition1 = c("motor_cortex", "visual_cortex"),
#'     # analogous, but for the second group in A/B testing.
#'     # note that you may alternatively this to `NA` to indicate
#'     # "everything except values in values_condition1"
#'     values_condition2 = c("prefrontal_cortex"),
#'     # a vector/array of 0 or more column names in dataset$samples
#'     # that should be used as additional regression variables
#'     colname_additional_variables = c("batch", "age")
#'   )
#'
#'   # Assume that the sample metadata table contains a column "group" with values
#'   # "A", "B", "C", "D". The following code will create the contrast "A vs B,C,D"
#'   # by setting the second set of values to `NA`
#'   dataset = add_contrast(
#'     dataset,
#'     colname_condition_variable = "group",
#'     values_condition1 = "A",
#'     values_condition2 = NA
#'   )
#'
#'   ## The "group" column in the sample table is used for group definitions in
#'   ## "all group" filtering. We here create two contrasts based on different
#'   ## regression variables in the sample tabel.
#'   dataset = remove_contrasts(dataset) # optionally, remove previously defined contrasts
#'   dataset = add_contrast(
#'     dataset,
#'     colname_condition_variable = "genotype",
#'     values_condition1 = "control",
#'     values_condition2 = "knockout",
#'     colname_additional_variables = "batch"
#'   )
#'   dataset = add_contrast(
#'     dataset,
#'     colname_condition_variable = "genotype",
#'     values_condition1 = "control",
#'     values_condition2 = c("knockout", "mutant1"),
#'     # note that we have the flexibility to add different regression variables per contrast
#'     colname_additional_variables = "batch"
#'   )
#'   # print an overview of all contrasts
#'   print_contrasts(dataset)
#'   # apply typical MS-DAP pipeline
#'   dataset = analysis_quickstart(
#'     dataset,
#'     filter_min_detect = 0, ## if DDA, we might not require minimum MS/MS counts
#'     filter_min_quant = 3,
#'     filter_fraction_detect = 0,
#'     filter_fraction_quant = 0.75,
#'     filter_by_contrast = FALSE,      ## we only want filtering across all groups
#'     filter_min_peptide_per_prot = 2, ## 2 peptides per protein
#'     dea_algorithm = "deqms",
#'     norm_algorithm = c("vwmb", "modebetween_protein"),
#'     output_qc_report = TRUE,
#'     output_abundance_tables = TRUE,
#'     output_dir = "C:/temp",         ## you may set this to NA to skip the QC report
#'     output_within_timestamped_subdirectory = TRUE
#'   )
#' }
#' @param dataset your dataset. Make sure you've already imported sample metadata
#' @param colname_condition_variable sample metadata column name that should be used for the experimental condition. Typically, this is the "group" column. Should be any of the values in `user_provided_metadata(dataset$samples)`
#' @param values_condition1 array of values from column `colname_condition_variable` that are the first group in the contrast. Note that
#' @param values_condition2 analogous to `values_condition1`, but for the second group. Note that you can set this to NA to indicate "everything except values in values_condition1"
#' @param colname_additional_variables optionally, sample metadata column names that should be used as additional regression variables (only the subset of `user_provided_metadata(dataset$samples)`, NOT including the value provided as parameter `colname_condition_variable`)
#' @seealso [print_contrasts()] to print an overview of defined contrasts, [remove_contrasts()] to remove all current contrasts (and respective filtering and DEA results)
#' @export
add_contrast = function(dataset, colname_condition_variable, values_condition1, values_condition2, colname_additional_variables = NULL) {
  # validate dataset
  if(!is.list(dataset) || !"samples" %in% names(dataset) || !is.data.frame(dataset$samples)) {
    append_log("Dataset does not contain a sample metadata table (dataset$samples)", type = "error")
  }
  if(!all(c("sample_id", "exclude") %in% colnames(dataset$samples)) || !all(dataset$samples$exclude %in% c(TRUE, FALSE))) {
    append_log("Expected column 'sample_id' and boolean column 'exclude' in dataset$samples", type = "error")
  }
  if("contrasts" %in% names(dataset) && !is.list(dataset$contrasts)) {
    append_log("dataset$contrasts must be a list. Overwriting this variable", type = "warning")
    dataset = remove_contrasts(dataset) # inits a new list
  }
  if(!"contrasts" %in% names(dataset)) {
    dataset = remove_contrasts(dataset) # inits a new list
  }

  # validate condition colname
  valid_sample_metadata_columns = user_provided_metadata(dataset$samples)
  if(length(colname_condition_variable) != 1 || is.na(colname_condition_variable) || colname_condition_variable == "" || !colname_condition_variable %in% valid_sample_metadata_columns) {
    append_log(paste0(
      'Invalid colname_condition_variable parameter: this should be a sample metadata column name that indicates the experimental condition, but provided value "', colname_condition_variable,
      '" does not exist in the sample metadata table! Note that matching is case-sensitive\nAvailable column names are: ',
      paste(valid_sample_metadata_columns, collapse=",")), type = "error")
  }
  tmp = unlist(dataset$samples[,colname_condition_variable], recursive = FALSE, use.names = FALSE)
  if(anyNA(tmp) || !is.character(tmp) || any(tmp == "")) {
    append_log('The sample metadata column provided in parameter colname_condition_variable must contain character values and not any NA or empty strings', type = "error")
  }

  # validate additional variables colnames
  if(!is.null(colname_additional_variables)) {
    if(anyNA(colname_additional_variables) || !is.character(colname_additional_variables) || any(colname_additional_variables %in% c("", "condition")) || anyDuplicated(colname_additional_variables)) {
      append_log('Invalid colname_additional_variables parameter: must be a character array. Cannot be NA, empty string or "condition". Duplicated not allowed.  (leave to default NULL if there are no additional regression variables)', type = "error")
    }
    if(colname_condition_variable %in% colname_additional_variables) {
      append_log('Parameter colname_additional_variables should not overlap with parameter colname_condition_variable', type = "error")
    }
    tmp = setdiff(colname_additional_variables, valid_sample_metadata_columns)
    if(length(tmp) > 0) {
      append_log(paste0(
        'Invalid colname_additional_variables parameter: this should be an array of sample metadata column names that indicate additional regression variables, but provided value(s) "',
        paste(tmp, collapse = ","),
        '" do not exists in the sample metadata table! Note that matching is case-sensitive\nAvailable column names are: ',
        paste(valid_sample_metadata_columns, collapse=",")), type = "error")
    }
  }

  # sorted sample metadata table only for samples in current contrast + transform condition into [1,2] classification
  samples = dataset$samples %>%
    filter(exclude == FALSE) %>%
    # importantly, select relevant regression variable columns **in order of priority** (with sample identifiers in front)
    select(sample_id, condition = !!as.symbol(colname_condition_variable), tidyselect::all_of(colname_additional_variables))

  # validate
  if(length(values_condition1) == 0 || anyNA(values_condition1) || !is.character(values_condition1) || any(values_condition1 == "" | !values_condition1 %in% unique(samples$condition)) || anyDuplicated(values_condition1)) {
    append_log(sprintf('parameter values_condition1 contains invalid values. Must match a value in sample metadata column "%s" (conform your provided column name specified in parameter colname_condition_variable). Must not be NA nor empty strings, no duplicates allowed', colname_condition_variable), type = "error")
  }
  if(length(values_condition2) == 1 && is.na(values_condition2)) { # infer second group
    values_condition2 = setdiff(unique(samples$condition), values_condition1)
    if(length(values_condition2) == 0) {
      append_log(paste0('parameter values_condition2 is set to NA but all unique values in sample metadata column "', colname_condition_variable, '" are already provided at parameter values_condition1 (note that "exclude" samples are disregarded)'), type = "error")
    }
  } else {
    if(length(values_condition2) == 0 || anyNA(values_condition2) || !is.character(values_condition2) || any(values_condition2 == "" | !values_condition2 %in% unique(samples$condition)) || anyDuplicated(values_condition2)) {
      append_log(sprintf('parameter values_condition2 contains invalid values. Must match a value in sample metadata column "%s" (conform your provided column name specified in parameter colname_condition_variable). Must not be NA nor empty strings, no duplicates allowed', colname_condition_variable), type = "error")
    }
  }
  tmp = intersect(values_condition1, values_condition2)
  if(length(tmp) > 0) {
    append_log(paste0('there must not be any overlap between parameters values_condition1 and values_condition2. Redundant values: ', paste(tmp, collapse = ",")), type = "error")
  }

  samples = samples %>%
    filter(condition %in% c(values_condition1, values_condition2)) %>%
    mutate(condition = as.integer((condition %in% values_condition2) + 1)) %>%
    arrange(condition)

  if(sum(samples$condition == 1L) < 2) {
    append_log(paste0('less than 2 samples (that are not set to "exclude") in sample metadata column "', colname_condition_variable, '" match the values_condition1 parameter; contrasts need at least 2 samples in each group'), type = "error")
  }
  if(sum(samples$condition == 2L) < 2) {
    append_log(paste0('less than 2 samples (that are not set to "exclude") in sample metadata column "', colname_condition_variable, '" match the values_condition2 parameter; contrasts need at least 2 samples in each group'), type = "error")
  }

  # convert each column to factor/numeric
  samples = enforce_sample_value_types(samples, redundant_columns = "error", show_log = TRUE)

  # check that each variable has condition-unique values
  for(col in colname_additional_variables) {
    y = samples %>% pull(!!col)
    if(n_distinct(y[samples$condition == 1L]) == 1 && n_distinct(y[samples$condition == 2L]) == 1) {
      append_log(paste0('"Invalid regression variable "', col, '" has no unique values in either condition  (solution: remove it from the colname_additional_variables parameter for the current contrast)'), type = "error")
    }
  }

  # compose label
  lbl_contrast = sprintf(
    "%s vs %s",
    paste(values_condition1, collapse=","),
    paste(values_condition2, collapse=",")
  )
  lbl = sprintf(
    "contrast: %s vs %s # condition_variable: %s",
    paste(values_condition1, collapse=","),
    paste(values_condition2, collapse=","),
    colname_condition_variable
  )
  if(length(colname_additional_variables) > 0) {
    lbl = paste0(lbl, " # additional_variables: ", paste(colname_additional_variables, collapse=","))
  }

  lbl = sub(" *$", "", lbl) # trim trailing whitespace
  append_log(lbl, type = "info")

  # double-check duplicates
  isdupe = FALSE
  for(i in seq_along(dataset$contrasts)) {
    if(dataset$contrasts[[i]]$label == lbl) {
      append_log(paste0("skipping a contrast that was already added to the dataset: ", lbl), type = "warning")
      isdupe = TRUE
      break
    }
  }

  if(!isdupe) {
    # finally, store the data
    dataset$contrasts[[length(dataset$contrasts) + 1]] = list(
      label = lbl,
      label_contrast = lbl_contrast,
      colname_condition_variable = colname_condition_variable,
      colname_additional_variables = colname_additional_variables,
      values_condition1 = values_condition1,
      values_condition2 = values_condition2,
      sampleid_condition1 = samples %>% filter(condition == 1) %>% pull(sample_id),
      sampleid_condition2 = samples %>% filter(condition == 2) %>% pull(sample_id),
      # the samples table is essential: we've added a 'condition' column and included additional variables relevant for the (user-specified) linear regression
      sample_table = samples,
      model_matrix = stats::model.matrix(~ . , data = samples %>% select(-sample_id)),
      # in the way the samples table is currently constructed,
      # we'd always need the "condition" regression variable as our main result
      regression_coefficient_name = "condition"
    )
  }

  return(dataset)
}



#' convert the data type of column in a sample metadata table to factors or numeric type
#'
#' @description importantly, columns "sample_index", "sample_id" and "shortname" are ignored alltogether
#' @param samples dataset$samples or a subset thereof. Typically, only columns sample_id and user-defined regression variables are included as columns
#' @param redundant_columns how to deal with columns in `samples` that are redundant with other columns. Options: "error" (default), "warning", "ignore"
#' @param show_log boolean, show data type per column in console? Default: `TRUE`
enforce_sample_value_types = function(samples, redundant_columns = "error", show_log = TRUE) {
  if(!is.data.frame(samples) || ncol(samples) == 0 || nrow(samples) == 0) {
    append_log("samples parameter must be a non-empty data.frame", type = "error")
  }
  if(length(redundant_columns) != 1 || !redundant_columns %in% c("error", "warning", "ignore")) {
    append_log('redundant_columns parameter must be any of these options: "error", "warning", "ignore"', type = "error")
  }
  if(length(show_log) != 1 || !show_log %in% c(TRUE, FALSE)) {
    append_log("show_log parameter must be either TRUE or FALSE", type = "error")
  }

  colnames_check = setdiff(colnames(samples), c("sample_index", "sample_id", "shortname"))
  for(col in colnames_check) {
    values = unlist(samples[,col], recursive = FALSE, use.names = FALSE)
    # input validation
    if(anyNA(values) || any(is.character(values) & values == "")) {
      append_log(paste0(
        "Invalid values (NA or empty string) in sample metadata (dataset$samples) column: ", col,
        "\n(solution: edit your sample metadata Excel table and replace empty values, then re-import it in R with import_sample_metadata() and finally re-run this function)"
      ), type = "error")
    }

    # transform data types to numeric/factor
    if(all(is.numeric(values) | is.factor(values))) {
      # nothing to do. i.e. we don't want to recast integers (which are valid is.numeric() nor factors)
      is_numeric = is.numeric(values)
    } else {
      values_recast = suppressWarnings(as.numeric(values))
      is_numeric = all(is.finite(values_recast))
      if(is_numeric) {
        samples[,col] = values_recast
      } else {
        values_recast = factor(values, levels = unique(values)) # don't resort, set levels according to input table sorting
        samples[,col] = values_recast
      }
    }

    # print data type
    if(show_log) {
      if(is_numeric) {
        append_log(paste0("numeric variable: ", col), type = "info")
      } else {
        append_log(paste0("categorical variable (R factor): ", col), type = "info")
      }
    }
  }

  # check that each variable is unique
  if(redundant_columns != "ignore") {
    # cast to numeric matrix so factors are checked as integer values (i.e. this helps spot patterns like A,A,B,B vs C,C,D,D - albeit not perfectly)
    m = as.matrix(
      samples %>%
        select(-tidyselect::any_of(c("sample_index", "sample_id", "shortname")) ) %>%
        mutate(across(tidyselect::where(is.factor), as.numeric))
    )
    # check if any of the previous columns is the same
    if(ncol(m) > 1) {
      for(i in 2L:ncol(m)) {
        for(j in 1L:(i-1L)) {
          if(all(m[,i] == m[,j])) {
            append_log(sprintf('Invalid regression variable "%s" is redundant, the same data is also encoded in regression variable "%s"  (solution: do not use this variable in your regression model/design)',
                               colnames(m)[i], colnames(m)[j]), type = redundant_columns)
          }
        }
      }
    }
  }

  return(samples)
}



#' Setup contrasts for downstream Differential Expression Analysis
#'
#' For more fine-grained control over specified contrasts, use the functions `remove_contrasts()` and `add_contrast()`.
#'
#' Note that a MS-DAP contrast for "A vs B" will return foldchanges for B/A in downstream output tables and data visualizations. For example, for the contrast "control vs disease" a positive log2 foldchange implies protein abundances are higher in the "disease" sample group.
#'
#' @param dataset your dataset. Make sure you've already imported sample metadata (so each sample is assigned to a sample group)
#' @param contrast_list a list that captures all contrasts that you want to compare. Check the examples for details.
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms in each statistical contrast tested downstream. Note that not all DEA algorithms may support this, consult documentation on individual methods for more info (start at `dea_algorithms()` )
#'
#' @examples
#' # a simple wild-type knockout study with only 2 groups, WT and KO
#' \dontrun{contrast_list = list(c("WT", "KO"))}
#' # multiple contrasts
#' \dontrun{contrast_list = list(c("cntrl", "group1"), c("cntrl", "group2"), c("group1", "group2"))}
#' # An example of multi-groups. The first contrast is just control-vs-group1
#' # as in the previous example. In the second contrast we compare control against
#' # multiple-groups combined as if they were one. Note the nested lists!
#' \dontrun{
#' contrast_list = list(c("control", "group1"), list("control", c("group1", "group2")))
#' }
#'
#' @export
setup_contrasts = function(dataset, contrast_list, random_variables = NULL) {
  if(!"samples" %in% names(dataset) || !is.data.frame(dataset$samples) || length(dataset$samples) == 0) {
    append_log("sample metadata table ('samples') is missing from the dataset. Run import_sample_metadata() prior to this function.", type = "error")
  }
  if(!"group" %in% colnames(dataset$samples)) {
    append_log("sample metadata table ('samples') must contain the column 'group'. Run import_sample_metadata() prior to this function.", type = "error")
  }
  if(length(contrast_list) == 0 || !is.list(contrast_list)) {
    append_log("contrast_list parameter must be a non-empty list", type = "error")
  }
  if(length(dataset$contrasts) > 0) {
    append_log("removing previous contrasts", type = "info")
  }
  dataset = remove_contrasts(dataset)

  random_variables = unique(random_variables)

  for(index in seq_along(contrast_list)) {
    contr = contrast_list[[index]]
    if(length(contr) != 2 || any(lengths(contr) == 0)) {
      append_log("each contrast should be a list with 2 elements (each is a non-empty array of sample group names)", type = "error")
    }

    dataset = add_contrast(dataset, colname_condition_variable = "group", values_condition1 = contr[[1]], values_condition2 = contr[[2]], colname_additional_variables = random_variables)
  }

  return(dataset)
}



#' add peptide and protein counts to samples tibble, so downstream code includes this in output tables and report figures relating to sample metadata
#'
#' @param peptides peptide tibble in long format
#' @param samples samples tibble in wide format
#' @param isdia boolean indicating this is a DIA dataset (DDA otherwise). Use is_dia_dataset(dataset) to obtain this from data
peptide_and_protein_counts_per_sample = function(peptides, samples, isdia) {
  peptide_counts = peptides %>%
    group_by(sample_id) %>%
    summarise(detected_peptides = sum(detect), all_peptides = n())

  protein_counts = peptides %>%
    group_by(sample_id, protein_id) %>%
    summarise(detect = any(detect)) %>%
    group_by(sample_id) %>%
    summarise(detected_proteins = sum(detect), all_proteins = n())

  i = match(samples$sample_id, peptide_counts$sample_id)
  j = match(samples$sample_id, protein_counts$sample_id)
  samples$detected_peptides = peptide_counts$detected_peptides[i]
  samples$detected_proteins = protein_counts$detected_proteins[j]

  if(!isdia) {
    samples$all_peptides = peptide_counts$all_peptides[i]
    samples$all_proteins = protein_counts$all_proteins[j]
  }

  return(samples)
}



#' Get column names that represent sample metadata provided by the user
#'
#' excludes all MS-DAP standard columns, like sample_id and shortname but also
#' computed metadata like peptide or protein counts
#'
#' @param samples e.g. dataset$samples
#' @export
user_provided_metadata = function(samples) {
  grep("^(sample_index|sample_id|shortname|exclude|contrast:.*|key_[a-z]+|((detected|quantified|all)_(peptides|proteins)))$", colnames(samples), value = T, invert = T, ignore.case = T)
}
