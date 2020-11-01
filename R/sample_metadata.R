
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

  dataset$samples = sample_metadata_from_file(sample_id = unique(dataset$peptides$sample_id), filename=filename)
  # add peptide and protein counts to samples tibble, so downstream code includes this in output tables and report figures relating to sample metadata
  dataset$samples = peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))
  return(dataset)
}



#' For advanced users: import sample metadata directly from the sample names using regular expressions
#'
#' @param dataset your dataset
#' @param sample_property_regex a named list of regular expressions to extract sample properties. name=sample property, value = array of length two with the regex and the replace string (first 2 arguments of sub())
#' @param sample_exclude_regex regex to select samples that should be flagged as 'exclude'
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

  s = gtools::mixedsort(unique(dataset$peptides$sample_id))
  openxlsx::write.xlsx(data.frame(sample_id=s, shortname=strip_common_substring(s), group="", exclude="FALSE", stringsAsFactors = F), file = filename)
}



#' placeholder title
#' @param sample_id todo
#' @param sample_property_regex todo
#' @param sample_exclude_regex todo
#' @param group_order todo
#' @param group_regex_array todo
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



#' placeholder title
#' @param sample_id todo
#' @param filename todo
#' @param group_order todo
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

  # read from disk
  if(is_type_xlsx) {
    df = openxlsx::read.xlsx(filename, sheet=1, sep.names="_")
  } else {
    df = as.data.frame(data.table::fread(filename))
  }
  if (! "sample_id" %in% colnames(df)) {
    append_log("sample_id is a required attribute for sample metadata", type = "error")
  }

  # input table should only contain samples matching this dataset
  samples_missing = setdiff(sample_id, df[,"sample_id"])
  samples_extra = setdiff(df[,"sample_id"], sample_id)
  if (length(samples_missing) > 0) {
    append_log(paste0("sample metadata table must contain all samples from this dataset. Missing samples:\n", paste(samples_missing, collapse="\n")), type = "error")
  }
  if (length(samples_extra) > 0) {
    append_log(paste("sample metadata table must contain only samples from this dataset. Additional samples that should be removed from metadata table:", paste(samples_extra, collapse=", ")), type = "error")
  }

  if(length(group_order) == 0 || any(!is.na(group_order))) {
    group_order = NA
  } else {
    group_order = unique(group_order)
  }

  # QC & sort
  return(sample_metadata_sort_and_filter(df, group_order = group_order))
}



#' placeholder title
#' @param df todo
#' @param sample_exclude_regex todo
#' @param group_order todo
#'
#' @importFrom gtools mixedorder
#' @importFrom stringr str_squish
sample_metadata_sort_and_filter = function(df, sample_exclude_regex = "", group_order = NA) {
  if (!is.data.frame(df) || length(df) == 0) {
    append_log("sample metadata must be a non-empty dataframe", type = "error")
  }

  ### first, format input table
  # convert all properties (/ column names) to lowercase
  colnames(df) = tolower(colnames(df))

  # dupes in column names
  if (anyDuplicated(colnames(df)) != 0) {
    print(colnames(df)[duplicated(colnames(df))])
    append_log("duplicated column names are not allowed", type = "error")
  }

  # replace whitespace in column names with underscores
  colnames(df) = gsub("\\s+", "_", colnames(df))

  # input validation on column names: valid characters
  cols_invalid = grepl("[^0-Z_ -.]", colnames(df), ignore.case = T)
  if (any(cols_invalid)) {
    append_log(paste("invalid characters in sample metadata column names; only alphanumeric characters and underscores are allowed. Invalid:", paste(colnames(df)[cols_invalid], collapse = ",")), type = "error")
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

  # replace NA values with empty string
  df[is.na(df)] = ""

  # trim whitespace, to prevent mistakes in input data (eg; group name has a trailing space in one row -->> should not be considered a distinct group)
  df = data.frame(apply(df, 2, stringr::str_squish), stringsAsFactors = F)

  # silently drop completely empty columns (eg; Excel may generate 'invisible' extra columns padded at the end)
  col_all_empty = colSums(df != "") == 0
  if(any(col_all_empty)) {
    # append_log(sprintf("removing empty columns from sample metadata table;", paste(colnames(df)[col_all_empty], collapse=", ")), type = "warning")
    df = df[, !col_all_empty, drop=F]
  }


  ### input validation on contents of the table
  # required properties/columns
  if (!all(c("sample_id", "shortname", "group") %in% colnames(df))) {
    append_log("sample_id, shortname and group are required attributes for sample metadata", type = "error")
  }

  # no empty strings, anywhere
  col_has_empty = colSums(df == "") > 0
  if (any(col_has_empty)) {
    append_log(paste("empty values are not allowed in sample metadata table. Columns with empty values;", paste(colnames(df)[col_has_empty], collapse=", ")), type = "error")
  }

  # dupes in sample_id
  if (anyDuplicated(df$sample_id) != 0) {
    print(df[duplicated(df$sample_id),])
    append_log("duplicated samples in the 'sample_id' column", type = "error")
  }


  if ("fractions" %in% colnames(df) && !"fraction" %in% colnames(df)) {
    append_log("a common typo in sample metadata column names is using 'fractions' instead of 'fraction', fixing now by renaming this column...", type = "warning")
    colnames(df)[colnames(df)=="fractions"] = "fraction"
  }

  # check for dupes in fraction * shortname
  if ("fraction" %in% colnames(df) && anyDuplicated(df[, c("shortname", "fraction")]) != 0) {
    print(df[duplicated(df[, c("shortname", "fraction")]),])
    append_log("duplicated fraction * sample-name ('shortname' column) combinations", type = "error")
  }

  if (!"fraction" %in% colnames(df) && anyDuplicated(df$shortname) != 0) {
    print(df[duplicated(df$shortname),])
    append_log("duplicated sample names in the 'shortname' column", type = "error")
  }

  # enforce valid characters in metadata
  rows_invalid = apply(df[, setdiff(colnames(df), "sample_id")], 1, function(x) any(grepl("[^0-Z_ -/\\.]", x, ignore.case = T))) # allow alphanumeric, underscore, space, minus
  if (any(rows_invalid)) {
    print(df[rows_invalid, , drop = F])
    append_log("invalid characters in sample metadata. outside of the column 'sample_id', only alphanumeric characters, underscores, - and spaces are allowed", type = "error")
  }
  rows_invalid = grepl("[^0-Z_+ -]", df$shortname, ignore.case = T)
  if (any(rows_invalid)) {
    print(df[rows_invalid, , drop = F])
    append_log("invalid characters in sample metadata; shortnames may only contain alphanumeric characters, underscores and spaces", type = "error")
  }


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
    # tolower is a string conversion (so works with both character and numeric input), %in% deals with NA's (although those are already caught upstream if all is well)
    df$exclude = tolower(df$exclude) %in% c("true", "1") # if user provides 0/1 flags
  }

  if (length(sample_exclude_regex) == 1 && !is.na(sample_exclude_regex) && sample_exclude_regex != "") {
    rows.remove = grepl(sample_exclude_regex, df$sample_id, ignore.case = T)
    if (any(rows.remove)) {
      # append_log(paste("samples on exclude list:", paste(df$sample_id[rows.remove], collapse=",")), type="info")
      df$exclude = rows.remove
    } else {
      append_log("your regex for removing unwanted samples did not match anything", type = "warning")
    }
  }

  ## backup from old systematics; groups must contain 2+ files
  # group_counts = table(df$group[df$exclude != TRUE])
  # fail_groups = names(group_counts)[group_counts < 2]
  # if (length(fail_groups) > 0) {
  #   print(df[df$group %in% fail_groups, , drop = F])
  #   append_log(paste("sample groups must have 2+ members (eg; a group of 1 has no replicates) -->> group(s): ", paste(fail_groups, collapse = ", ")), type = "error")
  # }

  # there must be at least 1 sample group with 1 non-exclude sample
  if(all(df$exclude == TRUE)) {
    append_log("there must be at least 1 sample that is not flagged as 'exclude'", type = "error")
  }

  # finally, add index as the first column
  df = data.frame(sample_index = 1:nrow(df), df, stringsAsFactors = F)

  return(as_tibble(df))
}


compose_contrast_name = function(a, b) {
  paste("contrast:", paste(unique(a), collapse = ","), "vs", paste(unique(b), collapse = ","))
}

decompose_contrast_name = function(x) {
  x = sub("^contrast: *", "", x)
  lapply(strsplit(x, split = " *vs *"), strsplit, split = " *, *")
}

#' Setup contrasts for downstream Differential Expression Analysis
#'
#' @param dataset your dataset, make sure you've already imported sample metadata (so each sample is assigned to a sample group)
#' @param contrast_list a list that captures all contrasts that you want to compare. Check the examples for details.
#' @param random_variables a vector of column names in your sample metadata table that are added as additional(!) regression terms in each statistical contrast tested downstream. Note that not all DEA algorithms may support this, consult documentation on individual methods for more info.
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
    append_log("sample metadata table ('samples') is missing from the dataset. Run import_sample_metadata() prior to this function.", type="error")
  }

  if(!"group" %in% colnames(dataset$samples)) {
    append_log("sample metadata table ('samples') must contain the column 'group'. Run import_sample_metadata() prior to this function.", type="error")
  }

  ## random variables
  random_variables = unique(random_variables)

  if(length(random_variables) > 0 && !all(random_variables %in% colnames(dataset$samples))) {
    append_log(paste("random_variables are case-sensitive and may only contain sample metadata (each value must match a column name in dataset$samples). Provided random_variables that are _not_ found in sample metadata table:", paste(setdiff(random_variables, colnames(dataset$samples)), collapse=", ")), type="error")
  }
  # not strictly needed, but help out the user by catching common mistakes; random effects as used in this pipeline must be 'additional metadata' besides the predefined contrasts
  ranvar_blacklist = c("sample_id", "group", "condition")
  if(any(ranvar_blacklist %in% random_variables)) {
    append_log(paste("random_variables should only contain 'additional metadata' that can be used to control for batch effects etc _besides_ the predefined contrasts in contrast_list, and _not_ contain these terms;", paste(intersect(random_variables, ranvar_blacklist), collapse=", ")), type="error")
  }

  # there can be no NA values or empty strings for sample metadata to be used as random variables in any downstream regression ('exclude' samples disregarded, these are never used in downstream DEA)
  if(length(random_variables) > 0) {
    ranvar_missing_values = NULL
    # iterate respective columns of sample metadata
    s = dataset$samples %>% filter(exclude == FALSE) %>% select(!!random_variables)
    for(v in colnames(s)) {
      x = s %>% pull(v)
      if(any(is.na(x) | x == "")) {
        # collect columns that contain missing values
        ranvar_missing_values = c(ranvar_missing_values, v)
      }
    }
    # report errors
    if(length(ranvar_missing_values) > 0) {
      append_log(paste("sample metadata table must not contain any missing values for columns that match the random_variables provided to this function. Sample metadata columns with missing values:", paste(ranvar_missing_values, collapse=", ")), type="error")
    }
  }

  # simply store list as dataset attribute
  # removes this variable (assign NULL) if no random_variables are supplied, as it should
  dataset$dea_random_variables = random_variables


  # reset cache
  dataset = invalidate_cache(dataset)

  column_contrasts = dataset_contrasts(dataset)
  # drop previous contrasts, if any
  tib = dataset$samples %>% select(-matches(column_contrasts))

  for (index in seq_along(contrast_list)) {
    contr = contrast_list[[index]]
    if (length(contr) != 2 || any(lengths(contr) == 0)) {
      append_log("each contrast should be a list with 2 elements (each is a non-empty array of sample group names)", type = "error")
    }

    # pretty-print labels
    groups_a = contr[[1]]
    groups_b = contr[[2]]
    lbl = compose_contrast_name(groups_a, groups_b)

    # input validation (is this a valid contrast list?)
    if (all(groups_a %in% groups_b) && all(groups_b %in% groups_a)) {
      append_log(paste("same group(s) on both sides of the", lbl), type = "error")
    }
    if (any(groups_a %in% groups_b) || any(groups_b %in% groups_a)) {
      append_log(paste("the same group cannot be on both sides of the", lbl), type = "error")
    }
    if (any(!groups_a %in% tib$group) || any(!groups_b %in% tib$group)) {
      append_log(paste("all groups in a contrast must be part of your dataset,", lbl), type = "error")
    }

    mask = as.integer(tib$group %in% groups_a)
    mask[tib$group %in% groups_b] = 2
    if ("exclude" %in% colnames(tib)) {
      mask[tib$exclude %in% TRUE] = 0
    }
    tib[, lbl] = as.integer(mask)

    ## QC: must have 2+ samples on both sides of the contrast
    count_samples_condition_a = sum(mask == 1)
    count_samples_condition_b = sum(mask == 2)
    if(count_samples_condition_a < 2 || count_samples_condition_b < 2) {
      append_log(sprintf("invalid contrast; '%s' vs '%s'. The former contains %d samples, the latter %d samples, while both should have at least 2 samples (disregarding samples that are flagged as 'exclude')",
                   paste(groups_a, collapse = ","), paste(groups_b, collapse = ","), count_samples_condition_a, count_samples_condition_b), type = "error")
    }

    append_log(lbl, type = "info")
  }
  dataset$samples = tib
  return(dataset)
}



#' add peptide and protein counts to samples tibble, so downstream code includes this in output tables and report figures relating to sample metadata
#' @param peptides todo
#' @param samples todo
#' @param isdia todo
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
