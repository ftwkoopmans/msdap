
#' Alternative FragPipe workflow, operates directly on a single psm.tsv file
#'
#' @description
#' For typical FragPipe workflows, `import_dataset_fragpipe_ionquant()` should be used instead of this function !
#'
#' This function directly parses data from a psm.tsv file and doesn't use MBR data produced by FragPipe;
#' use this function for legacy FragPipe data  OR  FragPipe runs where no experiment info was provided at all in
#' the FragPipe workflow tab.
#'
#' The most abundant PSM per peptide*sample is picked to represent abundance.
#'
#' @param filename the full file path of the input file
#' @param intensity_sum boolean value whether to take the sum intensity of all PSM per precursor (TRUE, default) or use the value from the PSM with highest abundance (FALSE)
#' @param acquisition_mode the type of experiment, should be a string. options: "dda" or "dia"
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param collapse_peptide_by if multiple data points are available for a peptide in a sample, at what level should these be combined? options: "sequence_modified" (recommended default), "sequence_plain", ""
#' @export
import_dataset_fragpipe_psm_file = function(filename, intensity_sum = TRUE, acquisition_mode, confidence_threshold = 0.01, collapse_peptide_by = "sequence_modified") {
  reset_log()
  append_log("importing FragPipe data from a single psm.tsv file using function import_dataset_fragpipe_psm_file(). If you only have 1 psm.tsv file from a FragPipe analysis without match-between-runs this is ok. But note that this function does not import IonQuant match-between-runs data!\nIn most cases you want to use MS-DAP function import_dataset_fragpipe_ionquant() instead which imports all IonQuant data including MBR abundance values (this MS-DAP function requires that in the FragPipe analysis, sample experiment info was provided. See further online documentation at MS-DAP GitHub > user guide > FragPipe)\n", type = "warning")

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

  check_parameter_is_numeric(confidence_threshold, minval_ = 0, maxval_ = 1)
  check_parameter_is_string(collapse_peptide_by)
  if(!(collapse_peptide_by %in% c("sequence_plain", "sequence_modified", ""))) {
    append_log('collapse_peptide_by parameter must be any of; "sequence_plain", "sequence_modified", ""', type = "error")
  }


  myfragpiperepair = function(DT) {
    # print(DT) # debug

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

    # replace non-finite intensities before working with intensity values downstream
    DT[!is.finite(intensity), intensity := 0]
    # sort such that the most abundant entry for each unique key comes out on top
    # also useful when we take the sum per group; e.g. rt will be value from highest intensity PSM
    data.table::setorder(DT, -intensity, na.last = TRUE)
    # set data.table grouping key
    data.table::setkey(DT, sequence_modified, charge, sample_id)

    # note that this parameter is from parent scope
    if(intensity_sum) {
      # replace all intensities per group by the within-group sum
      DT = DT[ , intensity := sum(intensity), by = .(sequence_modified, charge, sample_id)]
    } else {
      # to summarize peptide_id*sample*id by taking the PSM with highest abundance, we don't have to do anything
      # table is already sorted by 'highest intensity PSM first', so downstream code will pick our intended group representatives
    }

    # assign group index based on current table sorting, then select first entry per group
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
                                     return_decoys = FALSE,
                                     do_plot = FALSE)

  ds$acquisition_mode = acquisition_mode

  # convert RT to minutes
  ds$peptides$rt = ds$peptides$rt / 60


  ###### update protein identifiers, from protein table include the "Indistinguishable Proteins"
  file_proteins = paste0(dirname(filename), "/protein.tsv")
  if(!file.exists(file_proteins)) { # fallback
    file_proteins = paste0(dirname(filename), "/combined_protein.tsv")
  }

  if(file.exists(file_proteins)) {
    append_log(sprintf('parsing protein table "%s" that was found in same folder as psm.tsv to extract "Indistinguishable Proteins" (if this is undesired, move this file to another folder)', basename(file_proteins)), type = "info")

    tib_prot = parse_fragpipe_proteins_from_combined_proteins(file_proteins)
    # map data back to main peptide table. Note that above we selected the `Protein ID` column as proteingroup definition (these are the 'uniprot short ID' / only the accession)
    i = data.table::chmatch(ds$peptides$protein_id, tib_prot$protein_id)
    if(anyNA(i)) {
      append_log(sprintf('input data consistency problem; not all protein IDs found in psm.tsv file are present in column "Protein" from %s', basename(file_proteins)), type = "error")
    }
    ds$peptides$protein_id = tib_prot$protein_id_full[i]
  } else {
    append_log('no protein table (protein.tsv or combined_protein.tsv) available in same dir as psm.tsv, so cannot retrieve "Indistinguishable Proteins" from FragPipe (ambiguous protein IDs for each proteinGroup). So each proteinGroup is represented only by the leading protein ID.', type = "warning")
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
