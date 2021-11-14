
#' Import a label-free proteomics dataset from FragPipe presented in MSstats CSV format
#'
#' To generate output files that we here require in MS-DAP, configure FragPipe as follows:
#' - assign Experiment IDs in the workflow tab (you may simply set these all to 1)
#' - enable IonQuant ("Quant (MS1)" tab, as of FragPipe v15)
#' - optionally, enable match-between-runs
#'
#' This function will merge data from various FragPipe output files to construct a MS-DAP dataset
#'
#' @param path the full file path to the fragpipe output directory
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_fragpipe_ionquant = function(path, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()
  append_log("reading IonQuant data from FragPipe output folder...", type = "info")

  stopifnot(acquisition_mode %in% c("dda","dia"))
  # input validation; locate expected input files
  file_msstats = paste0(path, "/MSstats.csv")
  file_mbr = paste0(path, "/mbr_ion.tsv") # optional
  file_proteins = paste0(path, "/combined_protein.tsv")
  files_psm = dir(path, pattern = "^psm.tsv$", recursive = T, ignore.case = T, full.names = T)

  if (!file.exists(file_msstats)) {
    append_log(sprintf('Cannot find file "%s" at provided path "%s".\nPerhaps you forgot to enter experiment info in the FragPipe "Workflow" tab? That would prevent generation of MSstats.csv (as tested in FragPipe v15)', basename(file_msstats), path), type = "error")
  }
  if (!file.exists(file_proteins)) {
    append_log(sprintf('Cannot find file "%s" at provided path "%s"', basename(file_proteins), path), type = "error")
  }
  if (length(files_psm) == 0) {
    append_log(sprintf('no "psm.tsv" files found at provided path "%s"', path), type = "error")
  }


  ###### 1) process MSstats table into a MS-DAP formatted peptide table
  DT = data.table::fread(file_msstats)
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
    headers = unlist(strsplit(readLines(file_mbr, n = 1), "[\t,]+"))
    col_missing = setdiff(c("probability","intensity","apex_rt","modified_peptide","charge","acceptor_run_name"), headers)
    if (length(col_missing) > 0) {
      append_log(paste('Expected column names in mbr_ion.tsv table: probability, intensity, apex_rt, modified_peptide, charge, acceptor_run_name.\nmissing:', paste(col_missing, collapse = ", ")), type = "error")
    }

    # just to double-check and be robust to future changes, select per sample*peptide the unique entry with highest confidence (currently these entries are unique so the `distinct` call is not strictly needed, but safety first and cba about minor performance loss)
    tib_mbr = as_tibble(data.table::fread(file_mbr)) %>%
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
    append_log(sprintf('no match-between-runs data found (this file is missing: "%s"")', file_mbr), type = "info")
  }



  ###### 3) retention time and identification confidence from psm.tsv files
  for(filename in files_psm) {
    append_log(paste("Extract peptide retention time and identification confidence from:", filename), type = "info")

    # efficiently read peptide table, simplified from generic code @ import_dataset_in_long_format()
    headers = unlist(strsplit(readLines(filename, n = 1), "[\t,]+"))
    map_required = map_headers(headers, list(sample_id = "Spectrum File",
                                             sequence_plain = "Peptide",
                                             # sequence_modified = "", # by default, init this column as plain sequences because 'modified peptide' column is useless (see comments downstream),
                                             mz = "Calibrated Observed M/Z",
                                             modifications = "Assigned Modifications",
                                             charge = "Charge",
                                             rt = "Retention",
                                             confidence = "PeptideProphet Probability"),
                               error_on_missing = T, allow_multiple = F)

    DT = data.table::fread(filename, select = as.integer(map_required), check.names = T, stringsAsFactors = F)
    colnames(DT) = names(map_required)

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
      tib_modification_instructions = fragpipe_modification_to_instruction(unique(modseq_update$modifications))
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
  tib_prot = parse_fragpipe_proteins(file_proteins)
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



fragpipe_modification_to_instruction = function(modifications) {
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



#' Alternative FragPipe workflow, operates directly on psm.tsv file
#'
#' For typical FragPipe workflows, import_dataset_fragpipe_ionquant() should be used
#' instead of this function.
#'
#' This function directly parses data from the psm.tsv and doesn't use MBR data
#' presented in other files. The most abundant PSM per peptide*sample is picked to
#' represent abundance.
#'
#' @param filename the full file path of the input file
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @param do_plot logical indicating whether to create QC plots that are shown in the downstream PDF report (if enabled)
#' @export
import_dataset_fragpipe_psm_file = function(filename, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified", return_decoys = FALSE, do_plot = TRUE) {
  reset_log()
  stopifnot(acquisition_mode %in% c("dda","dia"))
  myfragpiperepair = function(DT) {
    # ### debug
    # print(head(DT))
    # hist(DT$confidence, breaks=50)
    # hist(log10(DT$confidence), breaks=50)
    # print(sum(DT$confidence >= 0.99))
    # print(sum(DT$confidence >= 0.99) / nrow(DT))
    # print(sum(DT$confidence >= 0.95))
    # print(sum(DT$confidence >= 0.95) / nrow(DT))

    ### specifically for fragpipe; strip full path + "interact-" + ".pep.xml" extension
    # example entry in input table; C:\DATA\fragpipe_test\interact-a05191.pep.xml
    DT[ , sample_id := gsub("(.*(\\\\|/)interact\\-)|(^interact\\-)|(\\.pep\\.xml$)", "", sample_id, ignore.case=F), by=sample_id]

    ### FragPipe output has "PeptideProphet Probability", which is 0~1 ranged score with 1 = 100% confident.
    # In this R package we use confidence score = qvalue (or pvalue). So as a band-aid solution we use 1-probability to approximate
    DT$confidence = 1 - DT$confidence

    ### FragPipe output lacks peptide sequences in modified sequence column if a peptide has no modifications
    rows = DT$sequence_modified == ""
    DT$sequence_modified[rows] = DT$sequence_plain[rows]
    DT$detect = is.finite(DT$confidence) & DT$confidence <= confidence_threshold

    ### summarize peptide_id*sample*id by taking the PSM with highest abundance
    data.table::setkey(DT, sequence_modified, charge, sample_id)
    # sort such that the most abundant entry for each unique key comes out on top
    data.table::setorder(DT, -intensity, na.last = TRUE)
    # simply select first entry
    DT = DT[ , index := seq_len(.N), by = .(sequence_modified, charge, sample_id)][index == 1]
    DT$index = NULL

    return(DT)
  }


  append_log("reading FragPipe PSM report...", type = "info")
  ds = import_dataset_in_long_format(filename,
                                     attributes_required = list(sample_id = "Spectrum.File",
                                                                protein_id = "Protein.ID",
                                                                sequence_modified = "Modified Peptide",
                                                                sequence_plain = "Peptide",
                                                                charge = "Charge",
                                                                rt = "Retention",
                                                                intensity = "Intensity",
                                                                confidence = "PeptideProphet.Probability"),
                                     select_unique_precursor_per_modseq = FALSE,
                                     remove_peptides_never_above_confidence_threshold = TRUE,
                                     confidence_threshold = confidence_threshold,
                                     custom_func_manipulate_DT_onload = myfragpiperepair,
                                     return_decoys = return_decoys,
                                     do_plot = do_plot)

  ds$acquisition_mode = acquisition_mode

  # convert RT to minutes
  ds$peptides$rt = ds$peptides$rt / 60


  ###### update protein identifiers, from protein table include the "Indistinguishable Proteins"
  file_proteins = paste0(dirname(filename), "/protein.tsv")
  if(!file.exists(file_proteins)) { # fallback
    file_proteins = paste0(dirname(filename), "/combined_protein.tsv")
  }

  if(file.exists(file_proteins)) {
    append_log(sprintf('parsing protein table "%s" to extract "Indistinguishable Proteins"', basename(file_proteins)), type = "info")

    tib_prot = parse_fragpipe_proteins(file_proteins)
    # map data back to main peptide table. Note that above we selected the `Protein ID` column as proteingroup definition. here, map to exact same column
    i = data.table::chmatch(ds$peptides$protein_id, tib_prot$protein_id_asis)
    if(any(is.na(i))) {
      append_log(sprintf('input data consistency problem; not all proteins in column "Protein ID" in %s are found in column "Protein" in %s', basename(filename), basename(file_proteins)), type = "error")
    }
    ds$peptides$protein_id = tib_prot$protein_id[i]
  } else {
    append_log('no protein table (protein.tsv or combined_protein.tsv) available in same dir as psm.tsv, cannot retrieve "Indistinguishable Proteins" (ambiguous protein IDs for each proteinGroup)', type = "warning")
    # if using the protein IDs directly from the `Protein ID` column in psm.tsv, enforce formatting of delimiters
    ds$peptides$protein_id = gsub("[, ;]+", ";", ds$peptides$protein_id)
  }


  ######

  # collapse peptides by plain or modified sequence (eg; peptide can be observed in some sample with and without modifications, at multiple charges, etc)
  if(collapse_peptide_by == "") {
    # if 'no collapse' is set, at least merge by modseq and charge. eg; there may be multiple peaks for some peptide with the same charge throughout retention time
    ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = "peptide_id")
    append_log("NOT collapsing peptides by plain/modified-sequence, thus 2 observations of the same sequence with different charge are considered a distinct peptide_id. This is not recommended for DDA!", type = "warning")
  } else {
    ds$peptides = peptides_collapse_by_sequence(ds$peptides, prop_peptide = collapse_peptide_by) # alternative, collapse modified sequences; prop_peptide = "sequence_modified"
  }

  log_peptide_tibble_pep_prot_counts(ds$peptides)
  return(ds)
}



#' parse a proteins.tsv or combined_protein.tsv from fragpipe output and return
#' a well formatted protein_id column that includes both "leading protein ID"
#' and "Indistinguishable Proteins"
#'
#' @param file_proteins full file path
parse_fragpipe_proteins = function(file_proteins) {
  tib_prot = as_tibble(data.table::fread(file_proteins))
  # input validation
  col_missing = setdiff(c("Protein","Protein ID","Indistinguishable Proteins"), colnames(tib_prot))
  if (length(col_missing) > 0) {
    append_log(sprintf('Expected column names in %s table: Protein, Protein ID, Indistinguishable Proteins.\nmissing:', file_proteins, paste(col_missing, collapse = ", ")), type = "error")
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
