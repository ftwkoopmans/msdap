
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
           isdecoy = FALSE, detect = FALSE, confidence = NA, rt = NA, is_mbr = FALSE)
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
      select(rt = apex_rt, sequence_modified = modified_peptide, charge, sample_id = acceptor_run_name) %>%
      distinct(sample_id, sequence_modified, charge, .keep_all = T) %>%
      mutate(id_sample_peptide = paste(sample_id, paste(sequence_modified, charge, sep="_")),
             rt = rt / 60) # convert RT to minutes

    # map back to main peptide table
    i = data.table::chmatch(tib_result$id_sample_peptide, tib_mbr$id_sample_peptide)
    tib_result$rt = tib_mbr$rt[i]
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
                                             sequence_modified = "Modified Peptide",
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

    tib_psm = as_tibble(DT) %>%
      mutate(sequence_modified_backup = sequence_modified, # for debugging / code review
             sequence_modified = sequence_plain, # start with plain sequence
             modification_list = strsplit(toupper(modifications), "[, ]+"),
             ### FragPipe output has "PeptideProphet Probability", which is 0~1 ranged score with 1 = 100% confident.
             # In this R package we use confidence score = qvalue (or pvalue). So as a band-aid solution we use 1-probability to approximate
             confidence = 1 - confidence,
             # convert RT to minutes
             rt = rt / 60,
             detect = is.finite(confidence) & confidence <= confidence_threshold
      )

    ## convert to a list of changes to be applied to plain sequences
    jobs = tibble(index = rep(1:nrow(tib_psm), lengths(tib_psm$modification_list)),
                  modification = unlist(tib_psm$modification_list)) %>%
      mutate(
        modification = gsub("(", "[", gsub(")", "]", modification, fixed = T), fixed = T), # replace brackets
        insert_string = gsub(".*(\\[.*)", "\\1", modification),
        site = gsub("^(\\d+).*", "\\1", modification)
      )

    # n-term and c-term modifications have a different notation, set their string position to begin/end of string
    rows = grepl("N-TERM", jobs$modification)
    jobs$site[rows] = 0
    jobs$insert_string[rows] = paste0("n", jobs$insert_string[rows])
    rows = grepl("C-TERM", jobs$modification)
    jobs$site[rows] = Inf
    jobs$insert_string[rows] = paste0("c", jobs$insert_string[rows])
    # as numerics
    jobs$site = as.numeric(jobs$site)
    # debug; print(jobs %>% distinct(modification, .keep_all=T), n=99)

    # make sure the largest string index/position is modified first (thereby retaining relative string positions)
    jobs = jobs %>% arrange(index, desc(site))

    ## apply string mutations
    # quite slow in R, non-optimized code but should be simple and robust. Optimize later when there are more users and test data available (to rule out edge-cases etc. for more efficient implementations)
    for(i in 1:nrow(jobs)) { # i=1
      # debug; tib_psm[jobs$index[i],]; jobs[i,]

      s = tib_psm$sequence_modified[ jobs$index[i] ]
      # modification all the way up front
      if(jobs$site[i] == 0) {
        s = paste0(jobs$insert_string[i], s)
      } else {
        # modification at the end
        if(!is.finite(jobs$site[i])) {
          s = paste0(s, jobs$insert_string[i])
        } else {
          # inject into sequence
          s = paste0(substring(s, 1, jobs$site[i]), jobs$insert_string[i], substring(s, jobs$site[i]+1))
        }
      }

      tib_psm$sequence_modified[ jobs$index[i] ] = s
    }

    ## map back to overall peptide tibble and transfer the observed retention time and identification confidence
    tib_psm$id_sample_peptide = paste(tib_psm$sample_id, paste(tib_psm$sequence_modified, tib_psm$charge, sep="_"))
    # i = index of an entry in current psm table @ overall peptide result table
    i = data.table::chmatch(tib_psm$id_sample_peptide, tib_result$id_sample_peptide)
    # don't overwrite values that we already pulled from the MBR table
    i[tib_result$is_mbr[i]] = NA
    # remove unmapped entries
    tib_psm = tib_psm[!is.na(i),]
    i = na.omit(i)
    # transfer data
    tib_result$rt[i] = tib_psm$rt
    tib_result$confidence[i] = tib_psm$confidence
    tib_result$detect[i] = tib_psm$detect
  }



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
