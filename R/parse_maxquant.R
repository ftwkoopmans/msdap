
#' Import a label-free proteomics dataset from MaxQuant
#'
#' @param path the directory that contains the search results (typically the 'txt' filter that contains files 'proteinGroups.txt' and 'evidence.txt')
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param remove_shared remove all shared peptides that could not be assigned to a single proteingroup as a unique- or razor-peptide ?
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
  file_evidence = path_exists(path, "evidence.txt", try_compressed = TRUE)

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
                             score_mbr = "Match score",
                             # mq_pep_id = "Peptide ID",
                             mq_mod_pep_id = "Mod. peptide ID")

  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  tibble_evidence = read_table_by_header_spec(file_evidence, attributes_required, attributes_optional, as_tibble_type = TRUE)


  rows_fail = !is.character(tibble_evidence$sequence_plain) | nchar(tibble_evidence$sequence_plain) <= 3
  if(any(rows_fail)) {
    append_log(sprintf("%d/%d rows (%.1f%%) in MaxQuant evidence.txt contain no plain sequence", sum(rows_fail), length(rows_fail), sum(rows_fail)/length(rows_fail) * 100), type = "error")
  }
  # bugfix: MaxQuant 2.1.1.0 may return empty values in "modified sequence" column sometimes (found some cases on MBR hits @ timsTOF DDA data)
  # We here try to repair using the "Mod. peptide ID" column
  rows_fail = !is.character(tibble_evidence$sequence_modified) | nchar(tibble_evidence$sequence_modified) <= 3
  if(any(rows_fail)) {
    if("mq_mod_pep_id" %in% colnames(tibble_evidence)) {
      # if peptide identifiers are available, borrow modified sequences from there
      if(all(is.numeric(tibble_evidence$mq_mod_pep_id) & is.finite(tibble_evidence$mq_mod_pep_id))) {
        # lookup table for the set of unique modified sequence IDs to modified sequence string
        tib_lookup = tibble_evidence %>%
          select(sequence_modified, mq_mod_pep_id) %>%
          distinct(mq_mod_pep_id, .keep_all = T) %>%
          mutate(fail = !is.character(sequence_modified) | sequence_modified == "")

        fraction_fail = sum(tib_lookup$fail) / nrow(tib_lookup)
        if(fraction_fail <= 0.1) {
          nrow_before = nrow(tibble_evidence)
          # replace all modified sequence using lookup table
          tibble_evidence = tibble_evidence %>%
            # join modified sequence from lookup table. Subset of valid sequences only
            select(-sequence_modified) %>%
            left_join(tib_lookup %>% filter(!fail), by = "mq_mod_pep_id") %>%
            # remove those entries in evidence table that failed to map
            filter(!is.na(sequence_modified))

          # print warning for data integrity issue and report on success rate of repairing
          append_log(sprintf("%d/%d rows (%.1f%%) in MaxQuant evidence.txt contain no modified sequence. After repairing %.2f%% of problematic peptides, we have to remove %d/%d peptides (%.2f%% of all rows in evidence.txt). Note that unrepairable entries are most likely MBR hits that were never identified by MS/MS, most likely this information is missing in evidence.txt due to a bug in MaxQuant (if input files for MS-DAP were MaxQuant output as-is)",
                             sum(rows_fail), length(rows_fail), sum(rows_fail)/length(rows_fail) * 100, (1-fraction_fail) * 100, sum(tib_lookup$fail), nrow(tib_lookup), (nrow_before - nrow(tibble_evidence)) / nrow_before * 100), type = "warning")
        } else {
          # halt if more than 10% of peptides would be removed
          append_log(sprintf("%d/%d rows (%.1f%%) in MaxQuant evidence.txt contain no modified sequence. We can only repair %.1f%% of peptides, so %d/%d peptides would be removed in total. Most likely a bug in MaxQuant, please report the mqpar.xml and evidence.txt to the MaxQuant team",
                             sum(rows_fail), length(rows_fail), sum(rows_fail)/length(rows_fail) * 100, (1-fraction_fail) * 100, sum(tib_lookup$fail), nrow(tib_lookup)), type = "error")
        }
      } else {
        # halt if we cannot repair due to incomplete modified peptide ID column
        append_log(sprintf("%d/%d rows (%.1f%%) in MaxQuant evidence.txt contain no modified sequence. Modified peptide ID column is available but contains invalid (non-numeric) values, so we cannot repair. Most likely a bug in MaxQuant, please report the mqpar.xml and evidence.txt to the MaxQuant team", sum(rows_fail), length(rows_fail), sum(rows_fail)/length(rows_fail) * 100), type = "error")
      }
    } else {
      # halt if we cannot repair due to lack of modified peptide ID column
      append_log(sprintf("%d/%d rows (%.1f%%) in MaxQuant evidence.txt contain no modified sequence. Modified peptide ID column also missing, so we cannot repair. Most likely a bug in MaxQuant, please report the mqpar.xml and evidence.txt to the MaxQuant team.", sum(rows_fail), length(rows_fail), sum(rows_fail)/length(rows_fail) * 100), type = "error")
    }
  }


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
#' @export
import_maxquant_proteingroups = function(path, remove_shared = T) {
  file_prot = path_exists(path, "proteinGroups.txt", try_compressed = TRUE)
  file_pep = path_exists(path, "peptides.txt", try_compressed = TRUE)

  ###### parse proteinGroups.txt
  ## read proteingroup file and validate input columns
  attributes_required = list(id = "id",
                             reverse = "Reverse",
                             majority_protein_ids = "Majority protein IDs",
                             peptide_id_is_razor = "Peptide is razor",
                             peptide_ids = "Peptide IDs")
  attributes_optional = list(contaminant = c("contaminant", "Potential contaminant"))

  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  MQ_prot = read_table_by_header_spec(file_prot, attributes_required, attributes_optional, as_tibble_type = TRUE)


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
    pepid = suppressWarnings(as.integer(unlist(l_pepid, use.names = F))),
    protein_id = rep(MQ_prot$majority_protein_ids, lengths(l_pepid)),
    is_razor = tolower(unlist(l_pepisrazor, use.names = F)) == "true"
  )
  if(any(!is.finite(pep2prot$pepid))) {
    append_log("Mapping between MaxQuant files, all peptide identifiers in the 'Peptide IDs' column in proteinGroups.txt must be integers", type = "error")
  }


  ###### collect peptide sequences for each peptide id
  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  MQ_pep = read_table_by_header_spec(file_pep, attributes_required = list(id = "id", sequence_plain = "Sequence"), as_tibble_type = TRUE)

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
#' @export
import_maxquant_peptides = function(file_peptides, remove_shared_peptides = T, pep2prot = NA) {
  reset_log()
  append_log("importing MaxQuant data from peptides.txt -->> we strongly recommend to import MaxQuant data using the function msdap::import_dataset_maxquant_evidencetxt() that used evidence.txt peptides.txt proteinGroups.txt", type = "warning")

  # will check for presence of file as well as .gz/.zip extension if file doesn't exist, will throw error if both do not exist
  file_peptides = path_exists(file_peptides, try_compressed = TRUE)

  ## read MaxQuant peptides.txt and validate input columns
  attributes_required = list(#id = "id",
                             reverse = "Reverse",
                             sequence_plain = "Sequence",
                             # confidence = "PEP",
                             protein_id = c("Leading razor protein", "Proteins"))
  attributes_optional = list(contaminant = c("contaminant", "Potential contaminant"),
                             proteingroups = c("Protein group IDs", "proteingroups"))

  # basically this reads the CSV/TSV table from file and maps column names to expected names.
  # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
  tibble_peptides = read_table_by_header_spec(file, attributes_required, attributes_optional, regex_headers = "^Intensity..+", as_tibble_type = TRUE)


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
