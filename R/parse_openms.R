

#' Import a label-free DDA proteomics dataset from OpenMS, experimental feature !
#'
#' This parser is based on output we generated with OpenMS 2.5 using input data and parameters/instructions from: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Tutorials/Data/iPRG2015_Hannes/
#'
#' Check the file "run_lfq.sh" in the above instructions (from OpenMS team) for exact workflow and settings used.
#'
#' By lack of more example mzTab files that we can adapt to, this code now **only** supports input data generated as in run_lfq.sh from the above example;
#'
#' The main dependency / complication is that we require the peptide- and protein-scores exactly as computed in the reference OpenMS workflow.
#' Thus, on each mzML file, apply search engine (eg; xtandem), PercolatorAdapter, FalseDiscoveryRate, IDFilter, IDScoreSwitcher.
#' finally run ProteomicsLFQ on the entire dataset.
#'
#' @param filename full path to the mzTab file
#'
#' @importFrom data.table fread setkey chmatch
#' @importFrom stringr str_sub
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @export
import_dataset_openms_mztab = function(filename) {
  # debug; test file; filename = "C:/DATA/iPRG2015_Hannes/iPRG2015_targeted_only.mzTab"
  # debug; test file; filename = "E:/DATA/PXD007683/openms/result.mzTab"
  reset_log()
  append_log("OpenMS data import is a new feature, please consult the documentations for limitations / work-in-progress", type = "warning")

  # first, check if input file exists
  check_parameter_is_string(filename)
  # will check for presence of file as well as .gz/.zip extension if file doesn't exist, will throw error if both do not exist
  filename = path_exists(filename, NULL, try_compressed = TRUE)


  decorate_modified_sequences = function(sequence, modifications) {
    result = sequence

    ## inject modifications by manual string manipulation
    # example: 3-UNIMOD:35,5-UNIMOD:35 = 2 modifications, to be inserted at positions 3 and 5 respectively (where 0 is start of the string, and 1 after the first character/amino-acid)
    i = which(!is.na(modifications) & !modifications %in% c("", "null", "na"))
    mods = strsplit(modifications[i], ",", fixed = T)
    for(index in seq_along(i)) { #index=381
      m = mods[[index]]
      mod_index = sub("-.*", "", m)

      pad = 0
      s = sequence[i[index]]
      for(z in seq_along(m)) { # z=1
        z_pos = pad + as.integer(mod_index[z])
        z_formatted = stringr::str_sub(m[z], nchar(mod_index[z])+2)
        z_size = nchar(z_formatted) + 2

        s = paste0(stringr::str_sub(s, 0, z_pos), "(", z_formatted, ")", stringr::str_sub(s, z_pos+1))
        # debug; rbind(s, paste0(stringr::str_sub(s, 0, z_pos), "(", z_formatted, ")", stringr::str_sub(s, z_pos+1)))

        pad = pad + z_size
      }
      result[i[index]] = s
    }

    return(result)
  }


  # read all lines as single column
  mztab_raw = data.table::as.data.table(data.frame(V1 = read_textfile_compressed(filename, as_table = F), stringsAsFactors = F) )

  # input validation: on the first 10 lines we should see the expected input filetype
  if(!all(c("MTD\tmzTab-version\t1.0.0", "MTD\tmzTab-type\tQuantification") %in% head(mztab_raw$V1, 10))) { # can we validate input mzTab type?
    append_log("OpenMS input file must be mzTab quantification format type, input requirements are strict, check the docs", type = "error")
  }

  # key = first 3 characters
  mztab_raw$headers = stringr::str_sub(mztab_raw$V1, 1, 3)
  data.table::setkey(mztab_raw, "headers")


  ### parse metadata
  raw_mtd = mztab_raw["MTD"][[1]]
  tib_mtd = tibble::as_tibble(data.table::fread(text = raw_mtd, header = F))
  mtd_is_percolator_prot_pep = any(grepl("protein_search_engine_score\\[1\\]\t.*(q-value|Percolator).*", raw_mtd, ignore.case = T)) &
    any(grepl("peptide_search_engine_score\\[1\\]\t.*(q-value|Percolator).*", raw_mtd, ignore.case = T))
  # example msTab line: MTD ms_run[11]-location  file:///home/sachsenb/OpenMS/openms-build/iPRG2015/JD_06232014_sample2_A.mzML
  filenames = tib_mtd %>% filter(grepl("ms_run\\[\\d+\\].{0,2}location", V2, ignore.case = T)) %>% pull(V3)
  # strip path and whitelisted extensions from filename
  filenames = gsub(regex_rawfile_strip_extension(), "", filenames, ignore.case=T)


  ###
  row_header_prot = unlist(strsplit(mztab_raw$V1[match("PRH", mztab_raw$headers)], "\t", fixed=T), use.names = F)
  row_header_pep = unlist(strsplit(mztab_raw$V1[match("PEH", mztab_raw$headers)], "\t", fixed=T), use.names = F)
  # row_header_psm = unlist(strsplit(mztab_raw$V1[match("PSH", mztab_raw$headers)], "\t", fixed=T), use.names = F)
  if(length(row_header_prot) < 2 || length(row_header_pep) < 2) { # can we find header lines?
    append_log("OpenMS mzTab file is missing protein or peptide headers (PRH PEH)", type = "error")
  }


  ### validate input
  if(length(filenames) < 2) {
    append_log("OpenMS mzTab file should have at least 2 filenames in MTD", type = "error")
  }
  if(any(duplicated(filenames))) {
    append_log("OpenMS mzTab file not have any duplicated filenames in MTD (we discard path, the filename must not be a dupe)", type = "error")
  }
  # only compatible with the OpenMS example for now
  # must have q-value or percolator scores for both peptide and protein. only support example workflow for now :/
  # below lines check required columns for proteins and peptides
  if(!mtd_is_percolator_prot_pep) {
    append_log("OpenMS mzTab file is missing percolator data. required are;  protein_search_engine_score\\[1\\]\t.*(q-value|Percolator)   and   peptide_search_engine_score\\[1\\]\t.*(q-value|Percolator)", type = "error")
  }
  if(!all(c("accession", "best_search_engine_score[1]") %in% row_header_prot)) {
    append_log("OpenMS mzTab file is missing columns from protein header, required are; accession, best_search_engine_score[1]", type = "error")
  }
  if(!all(c("accession", "sequence", "charge", "modifications", "search_engine_score[1]_ms_run[2]", "peptide_abundance_study_variable[2]", "opt_global_retention_time_study_variable[2]") %in% row_header_pep)) {
    append_log("OpenMS mzTab file is missing columns from peptide header, required are; accession, sequence, charge, modifications, search_engine_score[1]_ms_run[2], peptide_abundance_study_variable[2], opt_global_retention_time_study_variable[2]", type = "error")
  }


  ## parse protein and PSM data tables
  # TODO: relatively slow, use optimized code from above
  tib_prot = tibble::as_tibble(data.table::fread(text = mztab_raw["PRT"][[1]], header = F))
  colnames(tib_prot) = row_header_prot
  tib_pep = tibble::as_tibble(data.table::fread(text = mztab_raw["PEP"][[1]], header = F))
  colnames(tib_pep) = row_header_pep



  ### process data tables
  ## PRT
  tib_prot_filtered = tib_prot %>%
    rename(pep = `best_search_engine_score[1]`) %>%
    mutate(isdecoy = F, pep = suppressWarnings(as.numeric(pep))) %>%
    filter(is.finite(pep) & pep <= 0.05) %>% # protein-level FDR cutoff
    select(accession)

  # remove peptides not matching valid proteins + add mock isdecoy column
  tib_pep = tib_pep %>% filter(accession %in% tib_prot_filtered$accession) %>% add_column(isdecoy = F)

  ## PEP
  # temp peptide_id, functional but not so pretty; plain sequence + modifications + charge
  # for efficiency, pretty-print modified sequences after initial filtering
  tib_pep$peptide_id = paste0(tib_pep$sequence, "_", tib_pep$modifications, "_", tib_pep$charge)

  ## !! remove duplicate peptide*protein entries
  # these are typically peptides found at such distinct retention times that the OpenMS feature aligners could not link them
  tmp = paste0(tib_pep$peptide_id, tib_pep$accession)
  rows_fail = which(duplicated(tmp))
  if(length(rows_fail) > 0) {
    tib_pep = tib_pep %>% filter(!peptide_id %in% tib_pep$peptide_id[rows_fail])
  }

  ## precursors (peptide_id) that are listed multiple times; retain data from first line in table and use subsequent entries only to collect all protein accessions (don't have to do unique() because of upstream filters)
  pep2prot = tib_pep %>% group_by(peptide_id) %>% summarise(accession = paste(accession, collapse=";"))
  # drop duplicates
  tib_pep = tib_pep %>% distinct(peptide_id, .keep_all = T)
  # merge collapsed accessions back in
  tib_pep$accession = pep2prot$accession[data.table::chmatch(tib_pep$peptide_id, pep2prot$peptide_id)]


  ### update modified sequences and peptide_id
  tib_pep$sequence_modified = decorate_modified_sequences(tib_pep$sequence, tib_pep$modifications)
  tib_pep$peptide_id = paste0(tib_pep$sequence_modified, "_", tib_pep$charge)

  ### from wide to long format; confidence, intenstiy, rt
  cols_confidence = grep("search_engine_score\\[1\\]_ms_run\\[\\d+\\]", colnames(tib_pep), value = T)
  tib_pep_confidence = tidyr::pivot_longer(tib_pep %>% select(all_of(c("peptide_id", "accession", "sequence", "sequence_modified", "charge", "isdecoy", cols_confidence))), cols = -c("peptide_id", "accession", "sequence", "sequence_modified", "charge", "isdecoy"), names_to = "sample_id", values_to = "confidence") %>%
    mutate(confidence = suppressWarnings(as.numeric(confidence)), sample_id = stringr::str_sub(sample_id, 31, -2)) %>%
    filter(is.finite(confidence))

  # peptides must pass detect threshold in at least 1 sample
  peptide_id_valid = tib_pep_confidence %>%
    filter(confidence <= 0.01) %>% # peptide-level FDR cutoff
    distinct(peptide_id) %>%
    pull()
  cat(length(peptide_id_valid), "valid peptide_id after filters\n")
  tib_pep_filtered = tib_pep %>% filter(peptide_id %in% peptide_id_valid)

  #
  cols_intensity = grep("peptide_abundance_study_variable\\[\\d+\\]", colnames(tib_pep_filtered), value = T)
  tib_pep_filtered_intensity = tidyr::pivot_longer(tib_pep_filtered %>% select(all_of(c("peptide_id", cols_intensity))), cols = -peptide_id, names_to = "sample_id", values_to = "intensity") %>%
    mutate(intensity = suppressWarnings(as.numeric(intensity)), sample_id = stringr::str_sub(sample_id, 34, -2)) %>%
    filter(is.finite(intensity) & intensity > 0)

  #
  cols_rt = grep("opt_global_retention_time_study_variable\\[\\d+\\]", colnames(tib_pep_filtered), value = T)
  tib_pep_filtered_rt = tidyr::pivot_longer(tib_pep_filtered %>% select(all_of(c("peptide_id", cols_rt))), cols = -peptide_id, names_to = "sample_id", values_to = "rt") %>%
    mutate(rt = suppressWarnings(as.numeric(rt)), sample_id = stringr::str_sub(sample_id, 42, -2)) %>%
    filter(is.finite(rt))

  tib_pep_result = tib_pep_confidence %>%
    rename(sequence_plain = sequence,
           protein_id = accession) %>%
    dplyr::left_join(tib_pep_filtered_intensity, by=c("peptide_id", "sample_id")) %>%
    dplyr::left_join(tib_pep_filtered_rt, by=c("peptide_id", "sample_id")) %>%
    filter(is.finite(intensity) & is.finite(rt)) %>%
    mutate(sample_id = filenames[as.integer(sample_id)])


  tib_pep_result$rt = tib_pep_result$rt / 60
  tib_pep_result$intensity = log2(tib_pep_result$intensity)
  tib_pep_result$intensity[!is.na(tib_pep_result$intensity) & tib_pep_result$intensity < 0] = 0 # note; we already removed zero intensity values when importing. here, we threshold extremely low values

  tib_pep_result$detect = tib_pep_result$confidence <= 0.01 # peptide-level FDR cutoff


  tib_pep_result = peptides_collapse_by_sequence(tib_pep_result, prop_peptide = "sequence_modified")

  log_peptide_tibble_pep_prot_counts(tib_pep_result)
  return(list(peptides=tibble_peptides_reorder(tib_pep_result), proteins=empty_protein_tibble(tib_pep_result), acquisition_mode = "dda"))

  #### some test plots
  # hist(tib_pep_result$rt)
  # hist(log10(tib_pep_result$confidence))
  # plot(density(tib_pep_result$intensity))
  # lines(density(tib_pep_result$intensity[tib_pep_result$confidence <= 0.01]), col=2)
  # table(tib_pep_result$confidence <= 0.01)
  # print(tail(tib_pep_result, 100), n=100)

  #### backup code for parsing PSM table
  # col_prot_isdecoy = na.omit(grep("MS.1002217", row_header_prot, ignore.case = T, value = T))[1]
  # col_prot_pep = na.omit(grep("MS.1001493", row_header_prot, ignore.case = T, value = T))[1]
  # col_psm_isdecoy = na.omit(grep("MS.1002217", row_header_psm, ignore.case = T, value = T))[1]
  # col_psm_pep = na.omit(grep("MS.1001493", row_header_psm, ignore.case = T, value = T))[1]
  # col_psm_intensity = na.omit(grep("MS.(1000042|1001141|1001844)", row_header_psm, ignore.case = T, value = T)[1])
  #
  #
  # ### validate input
  # # must have at least 2 filenames
  # stopifnot(length(filenames) > 1)
  # stopifnot(!duplicated(filenames))
  # # check required columns for proteins and PSM
  # stopifnot(length(col_prot_isdecoy) == 1 && length(col_prot_pep) == 1 && "accession" %in% row_header_prot)
  # stopifnot(length(col_psm_isdecoy) == 1 && length(col_psm_pep) == 1 && length(col_psm_intensity) == 1 && all(c("accession", "PSM_ID", "sequence", "retention_time", "charge", "modifications") %in% row_header_psm) )
  #
  #
  # tib_psm = tibble::as_tibble(data.table::fread(text = mztab_raw["PSM"][[1]], header = F))
  # colnames(tib_psm) = row_header_psm
  #
  #
  # # PSM
  # tib_psm_filtered = tib_psm %>%
  #   rename(isdecoy = !!col_psm_isdecoy, pep = !!col_psm_pep, intensity = !!col_psm_intensity) %>%
  #   mutate(isdecoy = isdecoy %in% c(1, "1", TRUE, "true", "TRUE", "True"), pep = suppressWarnings(as.numeric(pep)), intensity = suppressWarnings(as.numeric(intensity)) ) %>%
  #   filter(accession %in% tib_prot_filtered$accession & !isdecoy & is.finite(intensity) & intensity > 0 & is.finite(pep) & pep <= 0.01) # PSM-level FDR cutoff
  #
  # tib_psm_filtered$run = as.integer(sub("ms_run\\[(\\d+)\\].*", "\\1", tib_psm_filtered$spectra_ref))
  # tib_psm_filtered$sample_id = filenames[tib_psm_filtered$run]
  #
  # ### end result tibble
  # tib_result = x %>% select(-PSM_ID, -modifications)
}
