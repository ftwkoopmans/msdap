
#' Test if the dataset contains DIA data, simply testing if the acquisition_mode equals "dia" (case insensitive)
#'
#' @param dataset dataset to test
#' @export
is_dia_dataset = function(dataset) {
  if(length(dataset$acquisition_mode) != 1) {
    append_log("'acquisition_mode' attribute missing from dataset", type = "error")
  }
  tolower(dataset$acquisition_mode) == "dia"
}


#' List the name of all contrasts in the samples table
#'
#' @param dataset a valid dataset
#' @export
dataset_contrasts = function(dataset) {
  if(!"samples" %in% names(dataset)) {
    append_log("invalid dataset, it lacks a samples table", type = "error")
  }
  grep("^contrast:", colnames(dataset$samples), ignore.case = T, value = T)
}


#' Print a short summary of a dataset to console
#'
#' @param dataset dataset to print
#' @export
print_dataset_summary = function(dataset) {
  ### report some basic stats
  cat(sprintf("%d samples (of which %d marked as 'exclude' by user) in %d sample groups; %s\n",
              nrow(dataset$samples), sum(dataset$samples$exclude), n_distinct(dataset$samples$group), paste(unique(dataset$samples$group), collapse=", ") ))

  if("intensity_all_group" %in% colnames(dataset$peptides)) {
    cat(sprintf("%d peptides from %d proteins, of which %d peptides pass user-defined filter criteria in all sample groups ('exclude' samples disregarded)\n",
                n_distinct(dataset$peptides$peptide_id), n_distinct(dataset$peptides$protein_id), n_distinct(dataset$peptides$peptide_id[!is.na(dataset$peptides$intensity_all_group)]) ))
  } else {
    cat(sprintf("%d peptides from %d proteins\n",
                n_distinct(dataset$peptides$peptide_id), n_distinct(dataset$peptides$protein_id) ))
  }

  column_contrasts = dataset_contrasts(dataset)

  ### present DEA results as a table
  if("de_proteins" %in% names(dataset) && is_tibble(dataset$de_proteins) && nrow(dataset$de_proteins) > 0 && length(dataset$proteins) > 0) {
    tib_dea = dataset$de_proteins %>%
      left_join(dataset$proteins, by="protein_id") %>%
      filter(!is.na(qvalue)) %>%
      arrange(qvalue) %>%
      group_by(contrast, dea_algorithm=algo_de) %>%
      summarise("#proteins tested" = n(),
                "signif @ user settings" = sum(signif),
                "qvalue <= 0.01" = sum(qvalue <= 0.01),
                "qvalue <= 0.05" = sum(qvalue <= 0.05),
                "top10 proteins at qvalue <= 0.05" = paste(head(gene_symbols_or_id[qvalue <= 0.05], 10), collapse=", ")
      ) %>%
      ungroup() %>%
      arrange(match(contrast, column_contrasts)) %>%
      mutate(contrast = sub("^contrast: *", "", contrast))
    cat("\ndifferential expression analysis, result summary:\n")
    print(remove_rownames(tib_dea), width=Inf, n=25) # print(as.data.frame(tib_dea), row.names=F, max = 25)
  }

  ### present DT results as a table
  if("dd_proteins" %in% names(dataset) && is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0 && length(dataset$proteins) > 0) {
    tib_report_diffdetects_summary = dataset$dd_proteins %>%
      filter(!is.na(zscore_count_detect)) %>%
      left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by="protein_id") %>%
      # sort such that top hits come first
      arrange(desc(abs(zscore_count_detect))) %>%
      # summary stats per contrast
      group_by(contrast) %>%
      summarise(`#proteins tested` = n(),
                `#abs(zscore) >= 3` = sum(abs(zscore_count_detect) >= 3),
                `top10` = tolower(paste(stringr::str_trunc(head(gene_symbols_or_id, 10), width = 10, side = "right"), collapse=", ") )) %>%
      ungroup() %>%
      arrange(match(contrast, column_contrasts)) %>%
      mutate(contrast = sub("^contrast: ", "", contrast))

    cat("\ndifferential detection analysis, result summary:\n")
    print(as.data.frame(tib_report_diffdetects_summary), row.names=F, max = 25)

    if(!is_dia_dataset(dataset)) {
      tib_report_diffdetects_summary = dataset$dd_proteins %>%
        filter(!is.na(zscore_count_quant)) %>%
        left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by="protein_id") %>%
        # sort such that top hits come first
        arrange(desc(abs(zscore_count_quant))) %>%
        # summary stats per contrast
        group_by(contrast) %>%
        summarise(`#proteins tested` = n(),
                  `#abs(zscore) >= 3` = sum(abs(zscore_count_quant) >= 3),
                  `top10` = tolower(paste(stringr::str_trunc(head(gene_symbols_or_id, 10), width = 10, side = "right"), collapse=", ") )) %>%
        ungroup() %>%
        arrange(match(contrast, column_contrasts)) %>%
        mutate(contrast = sub("^contrast: ", "", contrast))

      cat("\ndifferential detection analysis, counting 'quantified peptides' (includes MBR), result summary:\n")
      print(as.data.frame(tib_report_diffdetects_summary), row.names=F, max = 25)

    }
  }

}



#' Test if the minimum set of required information is present in a dataset
#'
#' @param dataset dataset to validate
#' @export
check_dataset_integrity = function(dataset) {
  property_required = c("peptides","proteins","samples", "acquisition_mode")
  if(!is.list(dataset) || !all(property_required %in% names(dataset))) {
    append_log(paste("dataset must be a list with at least the following properties:", paste(property_required, collapse=", ")), type = "error")
  }

  check_valid_tibble_peptides(dataset$peptides)
  check_valid_tibble_proteins(dataset$proteins)
  check_valid_tibble_samples(dataset$samples)
}



#' data integrity checks for peptide tibble
#'
#' @param tib peptide tibble (eg; dataset$peptides)
#' @export
check_valid_tibble_peptides = function(tib) {
  if(length(tib) < 2 || !is.data.frame(tib) || nrow(tib) == 0) {
    append_log("'peptides' variable must be a non-empty tibble", type = "error")
  }

  colname_missing = setdiff(c("sample_id", "protein_id", "peptide_id", "sequence_plain", "sequence_modified", "intensity", "detect", "isdecoy"), colnames(tib))
  if(length(colname_missing) > 0) {
    append_log(paste("peptides tibble is lacking required columns:", paste(colname_missing, collapse=", ")), type = "error")
  }

  if(any(!is.finite(tib$detect) | !is.logical(tib$detect))) {
    append_log("peptides tibble column 'detect' must only contain TRUE or FALSE", type = "error")
  }
  if(any(!is.finite(tib$isdecoy) | !is.logical(tib$isdecoy))) {
    append_log("peptides tibble column 'isdecoy' must only contain TRUE or FALSE", type = "error")
  }

  # character columns; test type and non-empty
  cols_character = c("sample_id", "protein_id", "peptide_id", "sequence_plain", "sequence_modified")
  colnames_fail = NULL
  for(col in cols_character) {
    x = tib %>% pull(!!col)
    if(any(is.na(x) | !is.character(x) | x == "")) {
      colnames_fail = c(colnames_fail, col)
    }
  }
  if(length(colnames_fail) > 0) {
    append_log(paste("these peptides tibble columns must be character type and non-empty:", paste(colnames_fail, collapse=", ")), type = "error")
  }

  # numeric columns; test type and non-infinite
  cols_numeric = grep("(^intensity)|(^rt$)", colnames(tib), ignore.case = T, value = T) # always at least 1, since intensity is required column
  colnames_fail = NULL
  for(col in cols_numeric) {
    x = tib %>% pull(!!col)
    if(!all(is.na(x) | (is.numeric(x) & is.finite(x)))) {
      colnames_fail = c(colnames_fail, col)
    }
  }
  if(length(colnames_fail) > 0) {
    append_log(paste("these peptides tibble columns must be numeric type and either NA or finite (eg; no Inf or -Inf):", paste(colnames_fail, collapse=", ")), type = "error")
  }

  return(TRUE)
  # # constraints; peptide_id*protein_id*sample_id must be unique
  # if(anyDuplicated(tib %>% select(peptide_id, protein_id, sample_id)) != 0) {
  #   append_log("peptides tibble; combinations of peptide_id*protein_id*sample_id must not contain any duplicates (eg; each peptide is in a sample just once)", type = "error")
  # }
}



#' data integrity checks for protein tibble
#'
#' @param tib protein tibble (eg; dataset$proteins)
#' @export
check_valid_tibble_proteins = function(tib) {
  if(length(tib) < 2 || !is.data.frame(tib) || nrow(tib) == 0) {
    append_log("'proteins' variable must be a non-empty tibble", type = "error")
  }

  colname_missing = setdiff(c("protein_id", "fasta_headers", "gene_symbols_or_id"), colnames(tib))
  if(length(colname_missing) > 0) {
    append_log(paste("proteins tibble is lacking required columns:", paste(colname_missing, collapse=", ")), type = "error")
  }

  # character columns; test type and non-empty
  cols_character = c("protein_id", "fasta_headers", "gene_symbols_or_id")
  colnames_fail = NULL
  for(col in cols_character) {
    x = tib %>% pull(!!col)
    if(any(is.na(x) | !is.character(x) | x == "")) {
      colnames_fail = c(colnames_fail, col)
    }
  }
  if(length(colnames_fail) > 0) {
    append_log(paste("these proteins tibble columns must be character type and non-empty:", paste(colnames_fail, collapse=", ")), type = "error")
  }

  # constraints; protein_id must be unique
  if(anyDuplicated(tib$protein_id) != 0) {
    append_log("proteins tibble column protein_id must not contain any duplicates", type = "error")
  }

  return(TRUE)
}



#' data integrity checks for sample tibble
#'
#' @param tib sample tibble (eg; dataset$samples)
#' @export
check_valid_tibble_samples = function(tib) {
  if(length(tib) < 2 || !is.data.frame(tib) || nrow(tib) == 0) {
    append_log("'samples' variable must be a non-empty tibble", type = "error")
  }


  colname_missing = setdiff(c("sample_id", "shortname", "group", "exclude"), colnames(tib))
  if(length(colname_missing) > 0) {
    append_log(paste("samples tibble is lacking required columns:", paste(colname_missing, collapse=", ")), type = "error")
  }

  if(any(!is.finite(tib$exclude) | !is.logical(tib$exclude))) {
    append_log("samples tibble column 'exclude' must only contain TRUE or FALSE", type = "error")
  }

  # character columns; test type and non-empty
  cols_character = c("sample_id", "shortname", "group")
  colnames_fail = NULL
  for(col in cols_character) {
    x = tib %>% pull(!!col)
    if(any(is.na(x) | !is.character(x) | x == "")) {
      colnames_fail = c(colnames_fail, col)
    }
  }
  if(length(colnames_fail) > 0) {
    append_log(paste("these samples tibble columns must be character type and non-empty:", paste(colnames_fail, collapse=", ")), type = "error")
  }

  # constraints
  if(anyDuplicated(tib$sample_id) != 0) {
    append_log("samples tibble column 'sample_id' must not contain any duplicates", type = "error")
  }
  ## below criterium only applies when there is no fractionation
  # if(anyDuplicated(tib$shortname) != 0) {
  #   append_log("samples tibble column 'shortname' must not contain any duplicates", type = "error")
  # }

  return(TRUE)
}



#' report number of detected precursor and plainseq, target and decoy
#' @param peptides peptide tibble in long format
log_peptide_tibble_pep_prot_counts = function(peptides) {
  report_unique = peptides %>% group_by(isdecoy) %>% summarise(n_precursor=n_distinct(peptide_id), n_seq=n_distinct(sequence_plain), n_prot=n_distinct(protein_id)) %>% arrange(isdecoy)
  append_log(sprintf("%d target precursors, %d (plain)sequences, %d proteins", report_unique$n_precursor[1], report_unique$n_seq[1], report_unique$n_prot[1]), type = "info")
  if(nrow(report_unique) > 1) { # if there are decoys
    append_log(sprintf("%d decoy precursors, %d (plain)sequences, %d proteins", report_unique$n_precursor[2], report_unique$n_seq[2], report_unique$n_prot[2]), type = "info")
  }
}



#' placeholder title
#' @param tib todo
tibble_peptides_reorder = function(tib) {
  tib %>% select(!!intersect(c("peptide_id", "protein_id", "sample_id", "sequence_plain", "sequence_modified", "confidence", "detect", "intensity", "rt"), colnames(tib)), everything())
}



#' placeholder title
#' @param peptides todo
empty_protein_tibble = function(peptides) {
  uprot = unique(peptides$protein_id)
  return(tibble(protein_id = uprot, fasta_headers = uprot, gene_symbols_or_id = uprot))
}
