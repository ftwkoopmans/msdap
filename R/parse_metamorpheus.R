
#' Import a label-free proteomics dataset from MetaMorpheus
#'
#' @param path the directory that contains the search results (eg; present files are AllProteinGroups.tsv, AllQuantifiedPeaks.tsv, etc.)
#' @param protein_qval_threshold qvalue threshold for accepting target proteins
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param only_unique_peptides remove all shared peptides? recommended to leave at FALSE
#' @export
import_dataset_metamorpheus = function(path, protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified", only_unique_peptides = FALSE) {
  append_log("reading MetaMorpheus search task...", type = "info")

  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }

  check_parameter_is_numeric(protein_qval_threshold)
  check_parameter_is_string(path)
  if (!dir.exists(path)) {
    append_log(paste("directory does not exist:", path), type = "error")
  }

  file_prot = path_exists(path, "AllProteinGroups.tsv")
  file_quantpeaks = path_exists(path, "AllQuantifiedPeaks.tsv")

  pep2prot = import_metamorpheus_proteingroups(file_prot, protein_qval_threshold = protein_qval_threshold, only_unique_peptides = only_unique_peptides)
  mm_peaks = import_metamorpheus_quantifiedpeaks(file_quantpeaks)
  # table(pep2prot$sequence_plain %in% mm_peaks$sequence_plain)
  # n_distinct(mm_peaks$full_sequence)
  # n_distinct(pep2prot$protein_id)


  # simply inner-join the peptide and protein tables
  # Note that upstream, we already enforced that each peptide_id is only assigned to one protein_id
  tib_result = mm_peaks %>%
    select(peptide_id, sample_id, sequence_modified, sequence_plain, charge, mz, rt, detect, intensity) %>%
    mutate(isdecoy = F) %>%
    inner_join(pep2prot %>% select(protein_id, sequence_plain), by = "sequence_plain")

  # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
  tib_result$intensity = log2(tib_result$intensity)

  # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }


  log_peptide_tibble_pep_prot_counts(tib_result)
  return(list(peptides=tibble_peptides_reorder(tib_result), proteins=empty_protein_tibble(tib_result), plots=list(), acquisition_mode = "dda"))
}



#' placeholder title
#' @param file_prot todo
#' @param protein_qval_threshold todo
#' @param only_unique_peptides todo
#'
#' @importFrom data.table fread
import_metamorpheus_proteingroups = function(file_prot, protein_qval_threshold = 0.05, only_unique_peptides = FALSE) {
  ############  read proteingroup file and validate input columns ############
  attributes_required = list(protein_id = "Protein.Accession",
                             qvalue = "Protein.QValue",
                             decoy_contaminant_target = "Protein.Decoy.Contaminant.Target",
                             count_peptides = "Number of Peptides",
                             unique_peptides = "Unique.Peptides",
                             shared_peptides = "Shared.Peptides")

  headers = unlist(strsplit(readLines(file_prot, n = 1), "\t"))
  col_indices = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)

  # workaround because MetaMorpheus output columns are bugged, trailing \t after header row :(
  mm_prot = as_tibble(data.table::fread(file_prot, select = as.numeric(col_indices), skip = 1, header = F, stringsAsFactors = F, quote = "", strip.white = T))
  colnames(mm_prot) = names(col_indices)

  # filter 'valid' results; target (not decoy or contaminant) + within FDR cutoff + has at least one unique peptide
  mm_prot = mm_prot %>%
    filter(is.finite(qvalue) & qvalue <= protein_qval_threshold & decoy_contaminant_target == "T") %>%
    arrange(desc(count_peptides))

  # MetaMorpheus seems to use a pipe symbol to delimit protein accessions, replace by semicolon for downstream compatability with our codebase
  mm_prot$protein_id = gsub("|", ";", mm_prot$protein_id, fixed = T)

  # bugfix for MetaMorpheus duplicated protein IDs, often with underscore appended numerical IDs. examples; G5EA30;Q92879;Q92879_2;Q92879_3;Q92879_4;Q92879   A0A0A0MRT6;B6VEX4;Q8IZP0;Q8IZP0_10;Q8IZP0_11;Q8IZP0_3;Q8IZP0_4;Q8IZP0_5;Q8IZP0_6;Q8IZP0_7;Q8IZP0_8;Q8IZP0
  # split all accessions by semicolon, then remove the appended numerical IDs (on the unlisted thing so we can use a vectorized regex. for better speed, use data.table), finally collapse the list by semicolon again
  l = strsplit(mm_prot$protein_id, ";", fixed = T)
  ul = sub("_\\d+$", "", unlist(l, recursive = F, use.names = F))
  l = relist(ul, l)
  mm_prot$protein_id = unlist(lapply(l, function(x) paste(unique(x), collapse = ";")), recursive = F, use.names = F)


  ############ peptide-to-protein mapping table ############
  l_unique = strsplit(mm_prot$unique_peptides, "|", fixed = T)
  l_shared = strsplit(mm_prot$shared_peptides, "|", fixed = T)
  names(l_unique) = names(l_shared) = mm_prot$protein_id
  # to tibble
  pep2prot = bind_rows(
    tibble(protein_id = rep(names(l_unique), lengths(l_unique)), sequence_plain = unlist(l_unique, use.names = F), is_unique = T),
    tibble(protein_id = rep(names(l_shared), lengths(l_shared)), sequence_plain = unlist(l_shared, use.names = F), is_unique = F)
  )

  if(only_unique_peptides) {
    cnt_prefilter = n_distinct(pep2prot$sequence_plain)
    pep2prot = pep2prot %>% filter(is_unique == T) %>% distinct(sequence_plain, .keep_all = T) # latter distinct() is not needed, but protects against upstream bugs
    cnt_postfilter = n_distinct(pep2prot$sequence_plain)
    append_log(sprintf("Parsing MetaMorpheus protein-groups, %d/%d peptides are not unique and therefore removed", cnt_prefilter-cnt_postfilter, cnt_prefilter), type = "info")
  } else {
    # because we sorted by largest protein-groups on top, we can here assign razor peptides simply by matching to first proteingroup
    pep2prot = pep2prot %>% distinct(sequence_plain, .keep_all = T)
    append_log("Parsing MetaMorpheus protein-groups, all non-unique peptides were assigned to their largest protein-group (by peptide count)", type = "info")
  }

  return(pep2prot)
}



#' placeholder title
#' @param file_quantpeaks todo
#'
#' @importFrom data.table fread
import_metamorpheus_quantifiedpeaks = function(file_quantpeaks) {
  attributes_required = list(sample_id = "File Name",
                             protein_id = "Protein Group",
                             full_sequence = "Full Sequence",
                             charge = "Precursor Charge",
                             rt = "Peak RT Apex",
                             mz = "Theoretical MZ",
                             detect_type = "Peak Detection Type",
                             intensity = "Peak intensity")
  headers = unlist(strsplit(readLines(file_quantpeaks, n = 1), "\t"))
  col_indices = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)

  # workaround because MetaMorpheus output columns are bugged, trailing \t after header row :(
  mm_peaks = as_tibble(data.table::fread(file_quantpeaks, select = as.numeric(col_indices), skip = 1, header = F, stringsAsFactors = F))
  colnames(mm_peaks) = names(col_indices)


  # force conversion of retention time and intensity to numeric -->> filter only quantified peaks -->> add detect type
  mm_peaks = mm_peaks %>%
    mutate(rt = suppressWarnings(as.numeric(rt)),
           mz = suppressWarnings(as.numeric(mz)),
           intensity = suppressWarnings(as.numeric(intensity)),
           detect = toupper(detect_type) == "MSMS") %>%
    filter(is.finite(rt) & rt != 0 & is.finite(intensity) & intensity != 0)

  # !! in some MetaMorpheus versions, the AllQuantifiedPeaks.tsv table contains multiple (modified) peptide sequences on one row in the "Full Sequence" column, delimited by a "|"
  # we just assume the first peptide is whatever has the highest PSM score and remove the rest
  mm_peaks %>%
    mutate(sequence_modified = gsub("\\|.*", "", full_sequence),
           sequence_plain = gsub("\\[.*?\\]", "", sequence_modified)) %>%
    unite(col = peptide_id, sequence_modified, charge, sep = "_", remove = FALSE)
           # peptide_id = paste(sequence_modified, charge, sep="_"))
}


### implementation v1; takes normalized peptide intensities from quantifiedpeptides table, respecting metamorpheus normalization. does not feature 'razor' peptide checks like v2 implementation, which is also simpler and more robust
# import_dataset_metamorpheus = function(path, protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified") {
#   append_log("reading MetaMorpheus search task...", type = "info")
#
#   if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
#     append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
#   }
#
#   check_parameter_is_numeric(protein_qval_threshold)
#   check_parameter_is_string(path)
#   if (!dir.exists(path)) {
#     append_log(paste("directory does not exist:", path), type = "error")
#   }
#
#   # check if all required files are present
#   file_prot = path_exists(path, "AllProteinGroups.tsv")
#   file_quantpeaks = path_exists(path, "AllQuantifiedPeaks.tsv")
#   # file_all_peptides = paste0(path, "/AllPeptides.psmtsv")
#   file_pep = paste0(path, "/AllQuantifiedPeptides_BaseSequences.tsv")
#   if (!file.exists(file_pep)) {
#     file_pep = paste0(path, "/AllQuantifiedPeptides.tsv")
#   }
#   if (!file.exists(file_pep)) {
#     append_log(sprintf('Cannot find file "AllQuantifiedPeptides_BaseSequences.tsv" or "AllQuantifiedPeptides.tsv" in directory "%s"', path), type = "error")
#   }
#
#
#   #### protein table
#   ## read proteingroup file and validate input columns
#   attributes_required = list(protein_id = "Protein.Accession",
#                              qvalue = "Protein.QValue",
#                              decoy_contaminant_target = "Protein.Decoy.Contaminant.Target",
#                              unique_peptides = "Unique.Peptides",
#                              shared_peptides = "Shared.Peptides")
#
#   headers = unlist(strsplit(readLines(file_prot, n = 1), "\t"))
#   col_indices = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
#
#   # only read columns of interest to speed up file parsing
#   # workaround because MetaMorpheus output columns are bugged, trailing \t after header row :(
#   mm_prot = as_tibble(data.table::fread(file_prot, select = as.numeric(col_indices), skip = 1, header = F, stringsAsFactors = F))
#   colnames(mm_prot) = names(col_indices)
#
#   # MetaMorpheus seems to use a pipe symbol to delimit protein accessions, replace by semicolon for downstream compatability with our codebase
#   mm_prot$protein_id = gsub("|", ";", mm_prot$protein_id, fixed = T)
#
#   # filter 'valid' results; target + qval cutoff
#   mm_prot = mm_prot %>% filter(is.finite(qvalue) & qvalue <= protein_qval_threshold & decoy_contaminant_target == "T")
#
#   ### peptide to protein
#   l_unique = strsplit(mm_prot$unique_peptides, "|", fixed = T)
#   l_shared = strsplit(mm_prot$shared_peptides, "|", fixed = T)
#   names(l_unique) = names(l_shared) = mm_prot$protein_id
#   # to tibble
#   pep2prot = bind_rows(
#     tibble(protein_id = rep(names(l_unique), lengths(l_unique)), sequence_plain = unlist(l_unique, use.names = F), is_unique = T),
#     tibble(protein_id = rep(names(l_shared), lengths(l_shared)), sequence_plain = unlist(l_shared, use.names = F), is_unique = F)
#   )
#   # fix protein id formats; remove underscore appended IDs from protein accessions. example; "P10636_1761"  @   C:\DATA\David\RAW\2019-03-20-06-03-53\Task2-SearchTask\AllQuantifiedPeaks.tsv
#   pep2prot$protein_id = gsub("_\\d*$", "", pep2prot$protein_id)
#   pep2prot = pep2prot %>% distinct(.keep_all = T)
#
#   ### remove shared peptides (any protein associated to multiple protein_id)
#   sequence_dupe = unique(pep2prot$sequence_plain[duplicated(pep2prot$sequence_plain)])
#   if (length(sequence_dupe) > 0) {
#     append_log(sprintf("%d/%d peptide plain sequences were shared among multiple proteinGroups, and consequentially removed", length(sequence_dupe), length(unique(pep2prot$sequence_plain))), type = "info")
#     pep2prot = pep2prot %>% filter(!sequence_plain %in% sequence_dupe)
#   }
#   # table(duplicated(pep2prot$sequence_plain))
#
#
#   #### the peptides file contains the peptide-level intensities, optionally normalized by MM
#   # workaround because MetaMorpheus output columns are bugged, trailing \t after header row :(
#   mm_pep = as_tibble(fread(file_pep, skip = 1, header = F), stringsAsFactors = F)
#   n = unlist(strsplit(readLines(file_pep, n = 1), "\t"))
#   colnames(mm_pep)[seq_along(n)] = make.names(n)
#   # validate input columns (throw error if some columns are missing)
#   map_headers(n, list(sequence_modified = "Sequence"), error_on_missing = T, allow_multiple = F)
#
#   # rename
#   colnames(mm_pep)[tolower(colnames(mm_pep)) == "sequence"] = "sequence_modified"
#   # duplicate into peptide_id as well
#   mm_pep$peptide_id = mm_pep$sequence_modified
#   # convert to plain sequence by regex
#   mm_pep$sequence_plain = gsub("\\[.*?\\]", "", mm_pep$sequence_modified)
#   # remove unused peptides
#   mm_pep = mm_pep %>% filter(sequence_plain %in% pep2prot$sequence_plain)
#
#   # extract detects in long format
#   tib_detect = mm_pep %>%
#     select(!!grep("(^peptide_id|Detection.Type..+)", colnames(mm_pep), ignore.case = T, value = T)) %>%
#     gather(sample_id, detect, -peptide_id, na.rm = T) %>%
#     mutate(sample_id = sub("^Detection.Type.", "", sample_id, ignore.case = T)) %>%
#     mutate(detect = grepl("^MSMS", detect, ignore.case = T))
#
#   # analogous for intensities, also collecting sequence_plain and sequence_modified
#   tib_pep = mm_pep %>%
#     select(!!grep("(^peptide_id|sequence_plain|sequence_modified|intensity..+)", colnames(mm_pep), ignore.case = T, value = T)) %>%
#     gather(sample_id, intensity, -peptide_id, -sequence_plain, -sequence_modified, na.rm = T) %>%
#     filter(is.finite(intensity) & intensity > 0) %>%
#     mutate(sample_id = sub("^intensity_", "", sample_id, ignore.case = T)) %>%
#     left_join(tib_detect, by = c("peptide_id", "sample_id")) %>%
#     mutate(detect = !is.na(detect) & detect == T)
#
#   # we can just fetch protein_id since we enfore a 1:n pep2prot mapping upstream
#   tib_pep$protein_id = pep2prot$protein_id[match(tib_pep$sequence_plain, pep2prot$sequence_plain)]
#   if(anyDuplicated(tib_pep[, c("sample_id", "peptide_id")]) != 0) {
#     append_log("sample_id,peptide_id pairs must be unique", type = "error")
#   }
#
#
#   #### get RT from quantified peptides file
#   attributes_required = list(sample_id = "File.Name",
#                              sequence_modified = "Full.Sequence",
#                              rt = "Peak.RT.Apex",
#                              intensity = "Peak.intensity")
#   headers = unlist(strsplit(readLines(file_quantpeaks, n = 1), "\t"))
#   col_indices = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
#
#   # only read columns of interest to speed up file parsing
#   # workaround because MetaMorpheus output columns are bugged, trailing \t after header row :(
#   mm_peaks = as_tibble(data.table::fread(file_quantpeaks, select = as.numeric(col_indices), skip = 1, header = F, stringsAsFactors = F))
#   colnames(mm_peaks) = names(col_indices)
#
#   # force conversion of retention time and intensity to numeric
#   mm_peaks$rt = suppressWarnings(as.numeric(mm_peaks$rt))
#   mm_peaks$intensity = suppressWarnings(as.numeric(mm_peaks$intensity))
#   # retention time foreach sequence_modified * sample_id
#   mm_peaks_rt = mm_peaks %>%
#     filter(is.finite(rt) & rt != 0 & is.finite(intensity) & intensity != 0) %>%
#     group_by(sample_id, sequence_modified) %>%
#     summarise(rt = mean(rt))
#
#   #### repair MetaMorpheus bug; filenames with a-typical characters, such as '-'are replaced with a . ONLY in the allpeptides file.  fix; clean filenames in both
#   tib_pep$sample_id = gsub("[^a-zA-Z0-9_]", "_", tib_pep$sample_id)
#   mm_peaks_rt$sample_id = gsub("[^a-zA-Z0-9_]", "_", mm_peaks_rt$sample_id)
#
#
#   tib_result = left_join(tib_pep, mm_peaks_rt, by = c("sample_id", "sequence_modified"))
#   # enforce detect type as boolean
#   tib_result$detect = !is.na(tib_result$detect) & tib_result$detect == T
#   # store intensities as log values (so we don't have to deal with integer64 downstream, which has performance issues)
#   tib_result$intensity = log2(tib_result$intensity)
#
#
#
#   if("isdecoy" %in% colnames(tib_result)) {
#     tib_result$isdecoy = !is.finite(tib_result$isdecoy) | tib_result$isdecoy == TRUE
#   } else {
#     tib_result$isdecoy = F
#   }
#   log_peptide_tibble_pep_prot_counts(tib_result)
#
#   # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
#   if(collapse_peptide_by == "") {
#     append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
#   } else {
#     tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
#   }
#   attr(tib_result, 'acquisition_mode') <- 'dda'
#
#   return(tib_result)
# }
