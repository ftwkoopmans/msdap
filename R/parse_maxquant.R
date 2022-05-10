
#' Import a label-free proteomics dataset from MaxQuant
#'
#' @param path the directory that contains the search results (typically the 'txt' filter that contains files 'proteinGroups.txt' and 'evidence.txt')
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param remove_shared remove all shared peptides that could not be assigned to a single proteingroup as a unique- or razor-peptide ?
#'
#' @importFrom data.table fread
#' @export
import_dataset_maxquant_evidencetxt = function(path, collapse_peptide_by = "sequence_modified", remove_shared = TRUE) {
  reset_log()
  append_log("reading MaxQuant 'txt' folder", type = "info")

  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }

  check_parameter_is_string(path)
  if (!dir.exists(path)) {
    append_log(paste("directory does not exist:", path), type = "error")
  }

  # first, check if input file exists
  file_evidence = path_exists(path, "evidence.txt")

  # get proteingroup information; accessions and mapping to peptides
  pep2prot = import_maxquant_proteingroups(path, remove_shared = remove_shared)

  ## read evidence file and validate input columns
  attributes_required = list(id = "id",
                             reverse = "Reverse",
                             sequence_plain = "Sequence",
                             intensity = "Intensity",
                             sequence_modified = "Modified sequence",
                             msms_count = "MS/MS Count",
                             calibrated_rt = "Calibrated retention time",
                             raw_file = "Raw file",
                             pep = "PEP",
                             charge = "Charge",
                             mz = "m/z")
  attributes_optional = list(contaminant = c("contaminant", "Potential contaminant"),
                             score_target = "Score",
                             score_mbr = "Match score")

  headers = unlist(strsplit(readLines(file_evidence, n = 1), "\t"))
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)

  # here we don't allow duplicate matches
  col_indices = c(map_required, map_optional)
  col_indice_dupe_names = names(col_indices)[col_indices %in% col_indices[duplicated(col_indices)]]
  if(length(col_indice_dupe_names) > 0) {
    append_log(paste('Duplicates encountered when finding column names in input data;', paste(col_indice_dupe_names, collapse = ", ")), type = "error")
  }

  # only read columns of interest to speed up file parsing
  tibble_evidence = as_tibble(data.table::fread(file_evidence, select = as.numeric(col_indices), check.names = T, stringsAsFactors = F))
  colnames(tibble_evidence) = names(col_indices)


  # remove contaminants
  col_contaminant = which(colnames(tibble_evidence) == "contaminant")
  if (length(col_contaminant) == 1) {
    tibble_evidence = tibble_evidence %>% filter(is.na(contaminant) | contaminant == "")
  } else {
    append_log("cannot find contaminant column in MaxQuant evidence.txt", type = "warning")
  }

  # remove reverse, peptides not in pep2prot mapping (eg; shared peptides/contaminant proteingroups, etc) and entries without quantification
  tibble_evidence = tibble_evidence %>%
    filter((is.na(reverse) | reverse == "") & sequence_plain %in% pep2prot$sequence_plain & is.finite(intensity) & intensity > 0) %>%
    unite(col = peptide_id, sequence_modified, charge, sep = "_", remove = FALSE)

  ## tibble_evidence$peptide_id = paste(tibble_evidence$sequence_modified, tibble_evidence$charge, sep = "_")
  # format essential data columns
  tibble_evidence$detect = is.finite(tibble_evidence$msms_count) & tibble_evidence$msms_count > 0
  tibble_evidence$calibrated_rt[!is.finite(tibble_evidence$calibrated_rt)] = NA
  tibble_evidence$pep[!is.finite(tibble_evidence$pep)] = NA
  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  tibble_evidence$intensity = log2(tibble_evidence$intensity)
  tibble_evidence$intensity[!is.na(tibble_evidence$intensity) & tibble_evidence$intensity < 0] = 0
  tibble_evidence$protein_id = pep2prot$protein_id[match(tibble_evidence$sequence_plain, pep2prot$sequence_plain)]
  # import raw confidence scores
  if(all(c("score_target", "score_mbr") %in% colnames(tibble_evidence)) && any(is.finite(tibble_evidence$score_target))) {
    tibble_evidence = tibble_evidence %>% mutate(confidence_score = score_target)
    rows = is.finite(tibble_evidence$score_mbr)
    if(any(rows)) {
      tibble_evidence$confidence_score[rows] = tibble_evidence$score_mbr[rows]
    }
    tibble_evidence$confidence_score[!is.finite(tibble_evidence$confidence_score)] = NA
  }

  if("isdecoy" %in% colnames(tibble_evidence)) {
    tibble_evidence$isdecoy = tibble_evidence$isdecoy %in% TRUE
  } else {
    tibble_evidence$isdecoy = F
  }

  # prep result tibble
  tib_result = tibble_evidence %>%
    select(peptide_id, protein_id, sample_id = raw_file, sequence_plain, sequence_modified, charge, mz, intensity, confidence = pep, isdecoy, rt = calibrated_rt, detect)
  if("confidence_score" %in% colnames(tibble_evidence)) {
    tib_result = tib_result %>% mutate(confidence_score = tibble_evidence$confidence_score)
  }

  # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  log_peptide_tibble_pep_prot_counts(tib_result)
  return(list(peptides=tibble_peptides_reorder(tib_result), proteins=empty_protein_tibble(tib_result), acquisition_mode = "dda"))
}



#' Import peptide-to-protein mappings from MaxQuant
#'
#' @param path the directory that contains the search results (typically the 'txt' filter that contains files 'proteinGroups.txt' and 'evidence.txt')
#' @param remove_shared remove all shared peptides that could not be assigned to a single proteingroup as a unique- or razor-peptide ?
#'
#' @importFrom data.table fread
#' @export
import_maxquant_proteingroups = function(path, remove_shared = T) {
  file_prot = path_exists(path, "proteinGroups.txt")
  file_pep = path_exists(path, "peptides.txt")

  ###### parse proteinGroups.txt
  ## read proteingroup file and validate input columns
  attributes_required = list(id = "id",
                             reverse = "Reverse",
                             majority_protein_ids = "Majority protein IDs",
                             peptide_id_is_razor = "Peptide is razor",
                             peptide_ids = "Peptide IDs")
  attributes_optional = list(contaminant = c("contaminant", "Potential contaminant"))

  headers = unlist(strsplit(readLines(file_prot, n = 1), "\t"))
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)

  # here we don't allow duplicate matches
  col_indices = c(map_required, map_optional)
  col_indice_dupe_names = names(col_indices)[col_indices %in% col_indices[duplicated(col_indices)]]
  if(length(col_indice_dupe_names) > 0) {
    append_log(paste('Duplicates encountered when finding column names in input data;', paste(col_indice_dupe_names, collapse = ", ")), type = "error")
  }

  # only read columns of interest to speed up file parsing
  MQ_prot = as_tibble(data.table::fread(file_prot, select = as.numeric(col_indices), check.names = T, stringsAsFactors = F))
  colnames(MQ_prot) = names(col_indices)


  # remove contaminants
  col_contaminant = which(colnames(MQ_prot) == "contaminant")
  if (length(col_contaminant) == 1) {
    MQ_prot = MQ_prot %>% filter((is.na(contaminant) | contaminant == "") & ((is.na(reverse) | reverse == "")))
  } else {
    append_log("cannot find contaminant column in MaxQuant proteinGroups.txt", type = "warning")
    MQ_prot = MQ_prot %>% filter((is.na(reverse) | reverse == ""))
  }


  # list; from protein_id to accessions
  ## extract peptide-id to proteingroup mapping, accounting for razor peptides
  l_pepid = strsplit(MQ_prot$peptide_ids, ";", fixed = T)
  l_pepisrazor = strsplit(MQ_prot$peptide_id_is_razor, ";", fixed = T)
  # must be of the exact same dimensions. eg; for every peptide id there must be a flag for 'is razor'
  if(!all(lengths(l_pepid) == lengths(l_pepisrazor))) {
    append_log("l_pepid and l_pepisrazor must be of the exact same dimensions", type = "error")
  }
  # to tibble
  pep2prot = tibble(
    pepid = unlist(l_pepid, use.names = F),
    protein_id = rep(MQ_prot$majority_protein_ids, lengths(l_pepid)),
    is_razor = tolower(unlist(l_pepisrazor, use.names = F)) == "true"
  )


  ###### collect peptide sequences for each peptide id
  # read peptide file
  MQ_pep = as.data.frame(fread(file_pep, select = c("id", "Sequence"), colClasses = "character"), stringsAsFactors = F)
  if(ncol(MQ_pep) != 2) {
    append_log("MaxQuant peptides.txt file must contain both columns 'id' and 'Sequence'", type = "error")
  }
  colnames(MQ_pep)[tolower(colnames(MQ_pep)) == "sequence"] = "sequence_plain"


  # align with pep2prot. must have a match for all peptide ids
  i = match(pep2prot$pepid, MQ_pep$id)
  if(any(is.na(i))) {
    append_log("Mapping between MaxQuant files, all peptide identifiers in the 'Peptide IDs' column in proteinGroups.txt must match an entry in the 'id' column of peptides.txt", type = "error")
  }
  pep2prot$sequence_plain = MQ_pep$sequence_plain[i]


  # remove shared peptides (not razor, and/or multiple protein_id associations) and create a n:1 peptide to protein mapping
  if (remove_shared) {
    pep2prot = pep2prot %>% group_by(sequence_plain)
    cnt_prefilter = dplyr::n_groups(pep2prot)

    pep2prot = pep2prot %>%
      filter(is_razor) %>%
      distinct(.keep_all = T) %>%
      filter(n() == 1)

    cnt_postfilter = dplyr::n_groups(pep2prot)
    pep2prot = pep2prot %>% ungroup()

    if(cnt_prefilter != cnt_postfilter) {
      append_log(sprintf("Parsing MaxQuant proteinGroups.txt, %d/%d peptides are not a razor peptide in any protein-group and therefore removed", cnt_prefilter-cnt_postfilter, cnt_prefilter), type = "info")
    }
  }

  return(pep2prot)
}



#' Import a label-free proteomics dataset from MaxQuant, using the peptides.txt file
#'
#' Recommended to import directly from the evidence file using \code{import_dataset_maxquant_evidencetxt} instead, as the peptides.txt file lacks retention-tine information !
#'
#' Provided mostly for compatability with incomplete datasets that lack the evidence.txt and proteinGroups.txt files
#'
#' @param file_peptides full file path to MaxQuant's peptides.txt output file
#' @param remove_shared_peptides remove all peptides that are assigned to multiple proteingroups (checked by having a semicolon in the 'Protein group IDs' column)
#' @param pep2prot if available, provide a 2 column table that provides a mapping from each peptide to a protein. Required columns; protein_id and sequence_plain. Set NA to use "Leading razor protein" as protein grouping
#'
#' @importFrom data.table fread
#' @export
import_maxquant_peptides = function(file_peptides, remove_shared_peptides = T, pep2prot = NA) {
  reset_log()
  append_log("importing MaxQuant data from peptides.txt -->> we strongly recommend to import MaxQuant data using the function msdap::import_dataset_maxquant_evidencetxt() that used evidence.txt peptides.txt proteinGroups.txt", type = "warning")

  ## read MaxQuant peptides.txt and validate input columns
  attributes_required = list(#id = "id",
                             reverse = "Reverse",
                             sequence_plain = "Sequence",
                             # confidence = "PEP",
                             protein_id = c("Leading razor protein", "Proteins"))
  attributes_optional = list(contaminant = c("contaminant", "Potential contaminant"),
                             proteingroups = c("Protein group IDs", "proteingroups"))

  headers = unlist(strsplit(readLines(file_peptides, n = 1), "\t"))
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)

  #
  map_intensity = grep("^Intensity..+", headers, ignore.case = T)
  names(map_intensity) = headers[map_intensity]
  if(length(map_intensity) == 0) {
    append_log("Cannot find any Intensity columns in MaxQuant peptides.txt input file", type = "error")
  }

  # here we don't allow duplicate matches
  col_indices = c(map_required, map_optional, map_intensity)
  col_indice_dupe_names = names(col_indices)[col_indices %in% col_indices[duplicated(col_indices)]]
  if(length(col_indice_dupe_names) > 0) {
    append_log(paste('Duplicates encountered when finding column names in input data;', paste(col_indice_dupe_names, collapse = ", ")), type = "error")
  }

  # only read columns of interest to speed up file parsing
  tibble_peptides = as_tibble(data.table::fread(file_peptides, select = as.numeric(col_indices), check.names = T, stringsAsFactors = F))
  colnames(tibble_peptides) = names(col_indices)

  if(!"confidence" %in% colnames(tibble_peptides)) {
    tibble_peptides$confidence = NA
  }


  ### filters
  # shared between protein-groups
  if("proteingroups" %in% colnames(tibble_peptides)) {
    if(remove_shared_peptides) {
      tibble_peptides = tibble_peptides[!grepl(";", tibble_peptides$proteingroups, fixed = T), ]
    }
    tibble_peptides = tibble_peptides %>% select(-proteingroups)
  }

  # remove contaminants
  col_contaminant = which(colnames(tibble_peptides) == "contaminant")
  if (length(col_contaminant) == 1) {
    tibble_peptides = tibble_peptides %>% filter(is.na(contaminant) | contaminant == "") %>% select(-contaminant)
  } else {
    append_log("cannot find contaminant column in MaxQuant peptides.txt", type = "warning")
  }
  # remove reverse
  tibble_peptides = tibble_peptides %>% filter((is.na(reverse) | reverse == "")) %>% select(-reverse)

  ## convert wide to long
  tib = tibble_peptides %>%
    pivot_longer(cols = c(-sequence_plain, -confidence, -protein_id), names_pattern = "^Intensity.(.+)", names_to = "sample_id", values_to = "intensity") %>%
    filter(intensity != 0) %>%
    add_column(isdecoy = F) #, -proteingroups
  tib$peptide_id = tib$sequence_modified = tib$sequence_plain


  ## remove peptides not in mapping table, then update protein IDs
  if(is_tibble(pep2prot) && nrow(pep2prot) > 0) {
    tib = tib %>% filter(sequence_plain %in% pep2prot$sequence_plain)
    tib$protein_id = pep2prot$protein_id[match(tib$sequence_plain, pep2prot$sequence_plain)]
  }


  ## conform to expected properties, standardized throughout input data parsers in this codebase
  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  tib$intensity = log2(tib$intensity)
  tib$intensity[!is.na(tib$intensity) & tib$intensity < 0] = 0 # note; we already removed zero intensity values when importing. here, we threshold extremely low values
  tib$rt = NA # no retention time info in peptides.txt
  tib$detect = T # cannot differentiate between MBR and MS/MS identified

  log_peptide_tibble_pep_prot_counts(tib)
  return(list(peptides=tibble_peptides_reorder(tib), proteins=empty_protein_tibble(tib), acquisition_mode = "dda"))
}
