
#' THIS FUNCTION SHOULD ONLY BE USED FOR LEGACY FRAGPIPE DATASETS THAT PRODUCED 'mbr_ion.tsv' OUTPUT FILES (e.g. FragPipe v15)
#'
#' To generate output files that we here require in MS-DAP, configure FragPipe as follows:
#' - assign Experiment IDs in the workflow tab (you may simply set these all to 1)
#' - enable IonQuant ("Quant (MS1)" tab, as of FragPipe v15)
#' - optionally, enable match-between-runs
#'
#' This function will merge data from various FragPipe output files to construct a MS-DAP dataset by
#' parsing PSM hits from "psm.tsv" files, collecting IonQuant peptide abundance values
#' from "MSstats.csv" and "mbr_ion.tsv", and finally fetching protein-group information from "combined_protein.tsv".
#'
#' @param path the full file path to the fragpipe output directory
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @importFrom stringi stri_sub_replace
#' @export
import_dataset_fragpipe_ionquant__legacy = function(path, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()
  append_log("using the legacy function for importing FragPipe data; import_dataset_fragpipe_ionquant__legacy(). If you have a legacy dataset from e.g. FragPipe v15, this is fine. But for more recent FragPipe datasets (e.g. that have 'combined_peptide.tsv' files in FragPipe output) you should use import_dataset_fragpipe_ionquant() instead", type = "warning")

  stopifnot(acquisition_mode %in% c("dda","dia"))
  # input validation; locate expected input files
  file_msstats = path_exists(path, "MSstats.csv", try_compressed = TRUE)
  file_mbr = path_exists(path, "mbr_ion.tsv", try_compressed = TRUE, silent = T) # optional
  file_proteins = path_exists(path, "combined_protein.tsv", try_compressed = TRUE)
  # check for nested files, regex analogous to path_exists()
  files_psm = dir(path, pattern = "^psm.tsv$", recursive = T, ignore.case = T, full.names = T)
  files_psm_compress1 = dir(path, pattern = "^psm.tsv.(zip|gz|bz2|xz|7z|zst|lz4)$", recursive = T, ignore.case = T, full.names = T)
  files_psm_compress2 = dir(path, pattern = "^psm.(zip|gz|bz2|xz|7z|zst|lz4)$", recursive = T, ignore.case = T, full.names = T)
  # de-duplicate using the path (so we don't have to deal with extensions in this de-dupe check)
  files_psm_dir = dirname(files_psm)
  files_psm_compress1_dir = dirname(files_psm_compress1)
  files_psm_compress2_dir = dirname(files_psm_compress2)
  files_psm = c(files_psm, files_psm_compress1[! files_psm_compress1_dir %in% files_psm_dir])
  files_psm = c(files_psm, files_psm_compress2[! files_psm_compress2_dir %in% files_psm_dir])

  # finally we can check if we found all required files
  if (file_msstats == "") {
    append_log(sprintf('Cannot find file "MSstats.csv" at provided path "%s".\nPerhaps you forgot to enter experiment info in the FragPipe "Workflow" tab? That would prevent generation of MSstats.csv (as tested in FragPipe v15)', path), type = "error")
  }
  if (file_proteins == "") {
    append_log(sprintf('Cannot find file "combined_protein.tsv" at provided path "%s"', path), type = "error")
  }
  if (length(files_psm) == 0) {
    append_log(sprintf('no "psm.tsv" files found at provided path "%s"', path), type = "error")
  }


  ###### 1) process MSstats table into a MS-DAP formatted peptide table
  DT = read_textfile_compressed(file_msstats, as_table = T)
  colnames(DT) = tolower(colnames(DT))
  # input validation
  col_missing = setdiff(c("proteinname","peptidesequence","precursorcharge","run","intensity"), colnames(DT))
  if (length(col_missing) > 0) {
    append_log(paste('Expected column names in MSstats table (not case sensitive): ProteinName, PeptideSequence, PrecursorCharge, ProductCharge, Run, Intensity.\nmissing:', paste(col_missing, collapse = ", ")), type = "error")
  }

  # efficient computation of plain sequences using data.table
  DT[ , sequence_plain := gsub("(\\[[^]]*\\])|(\\([^)]*\\))", "", peptidesequence), by=peptidesequence]

  # transform MSstats table into our standard format peptide tibble
  tib_result = as_tibble(DT) %>% select(sample_id = run, protein_id = proteinname, sequence_modified = peptidesequence, sequence_plain, charge = precursorcharge, intensity) %>%
    filter(is.finite(intensity) & intensity > 1) %>%
    mutate(intensity = log2(intensity),
           peptide_id = paste(sequence_modified, charge, sep="_"),
           id_sample_peptide = paste(sample_id, peptide_id),
           isdecoy = FALSE, detect = FALSE, confidence = NA, rt = NA, mz=NA, is_mbr = FALSE)
  rm(DT)



  ###### 2) MBR retention times from mbr_ion.tsv  (if available)
  if(file.exists(file_mbr)) {
    append_log(paste("Extract peptide retention time and identification confidence from:", file_mbr), type = "info")
    # input validation
    headers = unlist(strsplit( read_textfile_compressed(file_mbr, as_table = F, nrow = 1), "[\t,]+"))
    col_missing = setdiff(c("probability","intensity","apex_rt","modified_peptide","charge","acceptor_run_name"), headers)
    if (length(col_missing) > 0) {
      append_log(paste('Expected column names in mbr_ion.tsv table: probability, intensity, apex_rt, modified_peptide, charge, acceptor_run_name.\nmissing:', paste(col_missing, collapse = ", ")), type = "error")
    }

    # just to double-check and be robust to future changes, select per sample*peptide the unique entry with highest confidence (currently these entries are unique so the `distinct` call is not strictly needed, but safety first and cba about minor performance loss)
    tib_mbr = as_tibble(read_textfile_compressed(file_mbr, as_table = T, data.table = F)) %>%
      filter(is.finite(intensity) & intensity > 1) %>%
      arrange(desc(probability), desc(intensity)) %>%
      select(mz, rt = apex_rt, sequence_modified = modified_peptide, charge, sample_id = acceptor_run_name) %>%
      distinct(sample_id, sequence_modified, charge, .keep_all = T) %>%
      mutate(id_sample_peptide = paste(sample_id, paste(sequence_modified, charge, sep="_")),
             rt = rt / 60) # convert RT to minutes

    # map back to main peptide table
    i = data.table::chmatch(tib_result$id_sample_peptide, tib_mbr$id_sample_peptide)
    tib_result$rt = tib_mbr$rt[i]
    tib_result$mz = tib_mbr$mz[i]
    tib_result$is_mbr = !is.na(i)
  } else {
    append_log(sprintf('no match-between-runs data found. If you configured MBR in FragPipe, be aware that this data was NOT imported into MS-DAP because the "mbr_ion.tsv" file could not be found at: "%s"', path), type = "warning")
  }



  ###### 3) retention time and identification confidence from psm.tsv files
  for(filename in files_psm) {
    append_log(paste("Extract peptide retention time and identification confidence from:", filename), type = "info")

    # efficiently read peptide table, simplified from generic code @ import_dataset_in_long_format()
    attributes_required = list(sample_id = "Spectrum File",
                               sequence_plain = "Peptide",
                               # sequence_modified = "", # by default, init this column as plain sequences because 'modified peptide' column is useless (see comments downstream),
                               mz = "Calibrated Observed M/Z",
                               modifications = "Assigned Modifications",
                               charge = "Charge",
                               rt = "Retention",
                               confidence = "PeptideProphet Probability")

    # basically this reads the CSV/TSV table from file and maps column names to expected names.
    # (complicated) downstream code handles compressed files, efficient parsing of only the requested columns, etc.
    DT = read_table_by_header_spec(filename, attributes_required, as_tibble_type = F)

    ### specifically for fragpipe; strip full path + "interact-" + ".pep.xml" extension
    # example entry in input table; C:\DATA\fragpipe_test\interact-a05191.pep.xml
    DT[ , sample_id := gsub("(.*(\\\\|/)interact\\-)|(^interact\\-)|(\\.pep\\.xml$)", "", sample_id, ignore.case=F), by=sample_id]

    # init modified sequence like plain sequence
    DT[ , sequence_modified := sequence_plain]

    # add key for plainseq+mods
    DT[ , ismod := modifications != ""]
    DT[ ismod==T , id_plainseq_modifications := .GRP, by=.(sequence_plain, modifications)]

    # report if there are any sample mis-matches
    sample_missing = setdiff(unique(DT$sample_id), tib_result$sample_id)
    if (length(sample_missing) > 0) {
      append_log(paste('These samples from PSM data table are not in main MSstats.csv table;', paste(sample_missing, collapse = ", ")), type = "warning")
    }



    ### repair modified sequences !
    # By default, the format is different between psm.tsv and both (MSstats.csv and mbr_ion.tsv)
    # Example
    # MSstats.tsv: HQGVM[15.9949]VGM[15.9949]GQK
    # psm.tsv: "Modified Peptide" = "HQGVM[147]VGM[147]GQK" and in another column "Assigned Modifications" = "5M(15.9949), 8M(15.9949)"
    #
    # more examples from MSstats.csv (denoting n-terminal mods and fixed modifications)
    # n[42.0106]AYHSFLVEPISC[57.0215]HAWNK              psm.tsv = n[43]AYHSFLVEPISCHAWNK   note; fixed modifications seem absent from psm.tsv but not MSstats.csv
    # KLM[15.9949]DLEC[57.0215]SRDGLM[15.9949]YEQYR
    #
    # print( unique(unlist(strsplit(DT$modifications, "[, ]+"))) )

    # 1) unique set of peptide sequence + modifications (the set of strings to update).
    # modified peptides will have multiple entries in psm.tsv (e.g. same peptide identified in multiple samples), so taking unique set is important for efficiency
    modseq_update = as_tibble(DT) %>%
      filter(ismod) %>% # skip unmodified
      select(id_plainseq_modifications, sequence_plain, sequence_modified, modifications) %>%
      distinct(sequence_plain, modifications, .keep_all = T)

    if(nrow(modseq_update) > 0) {
      # 2) from a modification in FragPipe to a standardized instruction for updating the plain sequence such that it matches modified sequences in MSstats.csv
      tib_modification_instructions = fragpipe_modification_to_instruction__legacy(unique(modseq_update$modifications))
      # debug; tib_modification_instructions %>% filter(grepl("N-", input, fixed = T))

      # 3) string manipulations; translate every sequence + set of N modifications to an ordered list of substrings, then join them all (e.g. collapse with paste)
      for(i in 1:nrow(modseq_update)) { # i=grep(", N-", modseq_update$modifications, fixed=T)[1]
        tib = tib_modification_instructions %>% filter(input == modseq_update$modifications[i])
        stringi::stri_sub_replace(str=modseq_update$sequence_modified[i], from=tib$site+1, to=tib$site, value=tib$insert)
        for(j in 1:nrow(tib)) {
          if(is.na(tib$site[j])) {
            modseq_update$sequence_modified[i] = paste0(modseq_update$sequence_modified[i], tib$insert[j])
          } else {
            modseq_update$sequence_modified[i] = stringi::stri_sub_replace(str=modseq_update$sequence_modified[i], from=tib$site[j]+1, to=tib$site[j], value=tib$insert[j])
          }
        }
      }

      # 4) join with input table by `id_plainseq_modifications` key
      DT$sequence_modified = modseq_update$sequence_modified[match(DT$id_plainseq_modifications, modseq_update$id_plainseq_modifications)] # have to map from DT because of n:1 mapping to unique modseq
      DT[ is.na(sequence_modified) , sequence_modified := sequence_plain] # set plain sequence for rows without a 'modifications' entry
    }


    # 5) map back to overall peptide tibble and transfer the observed retention time and identification confidence
    DT[ , id_sample_peptide := paste(sample_id, paste(sequence_modified, charge, sep="_")) ]
    i = data.table::chmatch(DT$id_sample_peptide, tib_result$id_sample_peptide)
    # remove unmapped entries
    DT = DT[!is.na(i)]
    i = na.omit(i)
    # transfer data
    tib_result$mz[i] = DT$mz
    tib_result$rt[i] = DT$rt / 60 # convert RT to minutes
    # FragPipe output has "PeptideProphet Probability", which is 0~1 ranged score with 1 = 100% confident.
    # In this R package we use confidence score = qvalue (or pvalue). So as a band-aid solution we use 1-probability to approximate
    tib_result$confidence[i] = 1 - DT$confidence
  }


  # define 'detect'
  tib_result$detect = is.finite(tib_result$confidence) & tib_result$confidence <= confidence_threshold

  ###### 4) update protein identifiers, from protein table include the "Indistinguishable Proteins"
  tib_prot = parse_fragpipe_proteins__legacy(file_proteins)
  # map data back to main peptide table
  i = data.table::chmatch(tib_result$protein_id, tib_prot$protein_asis)
  if(any(is.na(i))) {
    append_log(sprintf('input data consistency problem; not all proteins in column "ProteinName" in %s are found in column "Protein" in %s', basename(file_msstats), basename(file_proteins)), type = "error")
  }
  tib_result$protein_id = tib_prot$protein_id[i]



  ###### 5) finally, construct MS-DAP dataset object

  ### collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  log_peptide_tibble_pep_prot_counts(tib_result)
  return(list(peptides=tibble_peptides_reorder(tib_result), proteins=empty_protein_tibble(tib_result), plots=list(), acquisition_mode = acquisition_mode))
}



#' THIS FUNCTION SHOULD ONLY BE USED FOR LEGACY FRAGPIPE DATASETS; HELPER FOR EXTRACTING PROTEIN MODIFICATIONS
#'
#' @param modifications characer vector
fragpipe_modification_to_instruction__legacy = function(modifications) {
  if(length(modifications) == 0) {
    return(tibble())
  }

  tibble(input = modifications) %>%
    mutate(mod = gsub("(", "[", gsub(")", "]", toupper(input), fixed = T), fixed = T), # replace brackets
           modlist = stringr::str_split(mod, "[, ]+")) %>%
    unnest(modlist) %>%
    mutate(
      # n-term and c-term modifications have a different notation, set their string position to begin/end of string
      modlist = sub("N-TERM", "0n", modlist, fixed = T),
      modlist = sub("C-TERM", "NAc", modlist, fixed = T),
      # now we can split the modification into a string that should be inserted and some position
      insert = sub("^(NA|\\d+)[A-Z]{0,1}([a-z]{0,1}\\[.*)$", "\\2", modlist),
      site = sub("^(NA|\\d+)[A-Z]{0,1}([a-z]{0,1}\\[.*)$", "\\1", modlist),
      site = suppressWarnings(as.numeric(site))
    ) %>%
    # make sure the largest string index/position is modified first (thereby retaining relative string positions)
    arrange(desc(site))
}



#' THIS FUNCTION SHOULD ONLY BE USED FOR LEGACY FRAGPIPE DATASETS; PARSER FOR PROTEIN TABLES
#'
#' parse a proteins.tsv or combined_protein.tsv from fragpipe output and return
#' a well formatted protein_id column that includes both "leading protein ID"
#' and "Indistinguishable Proteins"
#'
#' @param file_proteins full file path
parse_fragpipe_proteins__legacy = function(file_proteins) {
  tib_prot = as_tibble(read_textfile_compressed(file_proteins, as_table = T))
  # input validation
  col_missing = setdiff(c("Protein","Protein ID","Indistinguishable Proteins"), colnames(tib_prot))
  if (length(col_missing) > 0) {
    append_log(sprintf('Expected column names in %s table: Protein, Protein ID, Indistinguishable Proteins. Missing: %s', file_proteins, paste(col_missing, collapse = ", ")), type = "error")
  }

  tib_prot %>%
    # rename columns (lowercase and no spaces) and subset columns
    select(protein_asis = Protein, protein_id_asis = `Protein ID`, protein_indistinguishable = `Indistinguishable Proteins`) %>%
    # split 1:n protein identifiers in the `Protein ID` column
    mutate(protein_list = strsplit(protein_id_asis, "[, ;]+")) %>%
    # split 0:n protein identifiers in the `Indistinguishable Proteins` column
    mutate(protein_indistinguishable_list = strsplit(protein_indistinguishable, "[, ;]+")) %>%
    # unpack/unnest everything to long-format
    unnest(cols = protein_list, keep_empty=T) %>%
    unnest(cols = protein_indistinguishable_list, keep_empty=T) %>%
    # now we can apply vectorized functions; enforce short protein IDs
    mutate(protein_list = fasta_id_short(protein_list),
           protein_indistinguishable_list = fasta_id_short(protein_indistinguishable_list)) %>%
    # collapse everything again at the level of original keys
    group_by(protein_asis, protein_id_asis) %>%
    # note; relatively slow
    summarise(protein_id = paste(na.omit(unique(c(protein_list, protein_indistinguishable_list))), sep=";", collapse=";")) %>%
    ungroup()
}
