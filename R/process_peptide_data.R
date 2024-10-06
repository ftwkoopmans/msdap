
#' placeholder title
#' @param peptides todo
#' @param prop_peptide todo
peptides_collapse_by_sequence = function(peptides, prop_peptide = "sequence_plain") {
  tib_result = tibble_peptides_reorder(as_tibble(aggregate_tibble_by_datatables(peptides, prop_peptide)))
  return(tib_result)
}



#' Remove all peptides that match some protein-level filters
#'
#' @examples
#' # remove all gene symbols that start with krt (eg; keratins)
#' \dontrun{gene_symbols = "krt"}
#' #remove keratins and IGGs using regular expression against uniprot fasta headers
#' #(particularly useful for IP experiments);
#' \dontrun{regular_expression = "ig \\S+ chain|keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)"}
#'
#' @param dataset the dataset to filter
#' @param remove_irt_peptides try to find the irt spike-in peptides in the fasta file header. looking for; |IRT| and IRT_KIT and "Biognosys iRT" (case insensitive). default:FALSE
#' @param regular_expression careful here, regular expressions are powerful but complex matching patterns. A 'regex' that is matched against the fasta header(s) of a protein(group). case insensitive
#' @param gene_symbols an array of gene symbols that are to be matched against the fasta header(s) of a protein(group). case insensitive
#' @param print_nchar_limit max number of characters for the fasta headers (of removed proteins) that are shown in the log
#' @export
remove_proteins_by_name = function(dataset, remove_irt_peptides = FALSE, regular_expression = "", gene_symbols = c(), print_nchar_limit = 150) {
  # TODO: input validation

  # if(remove_not_present_in_fasta) {
  #   dataset$proteins = dataset$proteins %>% filter(protein_id %in% dataset$peptides$protein_id) # sync
  #   nprot = nrow(dataset$proteins)
  #   dataset$proteins = dataset$proteins %>% filter(!is.na(fasta_headers) & fasta_headers != protein_id) # actual filter
  #   nprot_post_filter = nrow(dataset$proteins)
  #   if(nprot_post_filter < nprot) {
  #     append_log(sprintf("%d / %d proteins had no mapping to fasta file and were thus removed", nprot - nprot_post_filter, nprot), type = "info")
  #   }
  # }

  regex_filter = NULL
  if(remove_irt_peptides) {
    regex_filter = c(regex_filter, "\\|IRT\\||IRT_KIT|Biognosys iRT")
  }

  if(length(print_nchar_limit) != 1 || !is.numeric(print_nchar_limit) || !is.finite(print_nchar_limit) || print_nchar_limit < 25) {
    append_log("print_nchar_limit parameter must be a single number (larger than 25)", type = "error")
  }
  print_nchar_limit = as.integer(ceiling(print_nchar_limit))

  if(!is.na(regular_expression) && regular_expression != "") {
    regex_filter = c(regex_filter, regular_expression)
  }

  gene_symbols = gene_symbols[!is.na(gene_symbols) & nchar(gene_symbols) > 1]
  if(length(gene_symbols) > 0) {
    if(any(grepl("[^[:alnum:]]", gene_symbols))) {
      append_log(paste("gene symbols should only contain alpha-numeric characters. gene_symbols:", paste(gene_symbols, collapse=" ")), type = "error")
    }
    regex_filter = c(regex_filter, sprintf(" GN=(%s)( |;|$)", paste(gene_symbols, collapse = "|")))
  }

  if(length(regex_filter) > 0) {
    # sync; subset protein tibble for those that are actually in the dataset
    dataset$proteins = dataset$proteins %>% filter(protein_id %in% dataset$peptides$protein_id)

    # wrap each regex in brackets, then combine using pipe symbol
    regex = paste(sapply(regex_filter, function(x) paste0("(",x,")")), collapse = "|")
    rows_remove = grepl(regex, dataset$proteins$fasta_headers, ignore.case = T)
    if(any(rows_remove)) {
      # store fasta headers that are about to be removed
      h = dataset$proteins$fasta_headers[rows_remove]
      # remove unwanted proteins from protein tibble, then sync peptide tibble
      dataset$proteins = dataset$proteins[!rows_remove, ]
      dataset$peptides = dataset$peptides %>% filter(isdecoy | protein_id %in% dataset$proteins$protein_id)
      # optionally, print the fasta headers for proteins that were removed
      charlim = print_nchar_limit
      h_too_long = nchar(h) > charlim
      h[h_too_long] = paste(substr(h[h_too_long], 1, charlim-4), "...")
      # status report & return results
      append_log(sprintf("%d/%d proteins matching your filters were removed from the dataset", sum(rows_remove), length(rows_remove)), type = "info")
      append_log(paste0("filtered proteins:\n  ", paste(h, collapse="\n  ")), type = "info")
      return(dataset)
    }
  }

  # default
  append_log("protein filters didn't match anything", type = "info")
  return(dataset)
}



#' Import fasta file(s) that match your dataset.
#'
#' For fasta files from uniprot, we extract the gene symbols for each protein-group (not available for non-uniprot fasta files).
#'
#' @param dataset your dataset
#' @param files an array of filenames, these should be the full path
#' @param fasta_id_type what type of fasta files are these? options: "uniprot" (highly recommended) or otherwise any other character string (as we have no special rules for generic fasta files)
#' @param protein_separation_character the separation character for protein identifiers in your dataset. Most commonly this is a semicolon (eg; in maxquant/metamorpheus/skyline/etc.)
#' @param uppercase_symbols convert all gene symbols to upper case? default: TRUE
#' @return table where
#' protein_id = provided proteingroup identifier,
#' accessions = result from fasta_id_short applied to each semicolon-delimited element in protein_id (result is a semicolon-collapsed string),
#' fasta_headers = analogous to accessions, but matching the full FASTA header to each 'accessions' element,
#' gene_symbols = full set of gene symbols, '-' where missing in FASTA, matching each element in 'accessions',
#' gene_symbols_or_id = unique set of valid 'gene_symbol', or the FASTA full/long ID when there is no gene information
#' gene_symbol_ucount = number of unique gene_symbols for this proteingroup (i.e. unique valid elements in 'gene_symbols')
#' @export
import_fasta = function(dataset, files = NULL, fasta_id_type = "uniprot", protein_separation_character = ";", uppercase_symbols = TRUE) {
  if(length(files) == 0) {
    append_log("no fasta files were provided, downstream analysis/code will use protein IDs as surrogate for fasta headers", type = "info")
    dataset$proteins = empty_protein_tibble(peptides)
  } else {
    dataset$proteins = import_protein_metadata_from_fasta(
      protein_id = dataset$peptides %>% filter({if("isdecoy" %in% names(.)) isdecoy else F} != TRUE) %>% distinct(protein_id) %>% pull(),
      fasta_files = files,
      fasta_id_type = fasta_id_type,
      protein_separation_character = protein_separation_character,
      uppercase_symbols = uppercase_symbols
    )
  }

  return(dataset)
}



#' Parse fasta files and construct a protein metadata table for provided protein identifiers
#'
#' @param protein_id an array of protein identifiers (that should be available in provided `fasta_files`)
#' @param fasta_files an array of filenames, these should be the full path
#' @param fasta_id_type what type of fasta files are these? options: "uniprot" (highly recommended) or otherwise any other character string (as we have no special rules for generic fasta files)
#' @param protein_separation_character the separation character for protein identifiers in your dataset. Most commonly this is a semicolon (eg; in maxquant/metamorpheus/skyline/etc.)
#' @param uppercase_symbols convert all gene symbols to upper case? default: TRUE
#' @return table where
#' protein_id = provided proteingroup identifier,
#' accessions = result from fasta_id_short applied to each semicolon-delimited element in protein_id (result is a semicolon-collapsed string),
#' fasta_headers = analogous to accessions, but matching the full FASTA header to each 'accessions' element,
#' gene_symbols = full set of gene symbols, '-' where missing in FASTA, matching each element in 'accessions',
#' gene_symbols_or_id = unique set of valid 'gene_symbol', or the FASTA full/long ID when there is no gene information
#' gene_symbol_ucount = number of unique gene_symbols for this proteingroup (i.e. unique valid elements in 'gene_symbols')
import_protein_metadata_from_fasta = function(protein_id, fasta_files, fasta_id_type = "uniprot", protein_separation_character = ";", uppercase_symbols = TRUE) {
  if(length(protein_id) == 0 || !all(!is.na(protein_id) & is.character(protein_id) & nchar(protein_id) > 0)) {
    append_log('protein_id parameter must be an array of non-empty strings (no NA values, no "")', type = "error")
  }
  # valid paths for `fasta_files` is checked downstream and includes guessing compressed files / filename variations
  if(length(fasta_files) == 0 || !all(!is.na(fasta_files) & is.character(fasta_files) & nchar(fasta_files) > 0)) {
    append_log('fasta_files parameter must be an array of non-empty strings (no NA values, no "")', type = "error")
  }
  if(length(fasta_id_type) != 1 || is.na(fasta_id_type) || !is.character(fasta_id_type) || nchar(fasta_id_type) == 0) {
    append_log('fasta_id_type parameter must be a non-empty string (no NA values, no "")', type = "error")
  }
  if(length(protein_separation_character) != 1 || is.na(protein_separation_character) || !is.character(protein_separation_character) || nchar(protein_separation_character) == 0) {
    append_log('protein_separation_character parameter must be a non-empty string (no NA values, no "")', type = "error")
  }
  if(!uppercase_symbols %in% c(TRUE, FALSE)) {
    append_log('uppercase_symbols parameter must be a boolean value (TRUE or FALSE)', type = "error")
  }


  protein_id = unique(protein_id)

  #### read data from fasta files into tibble
  fasta = fasta_parse_file(fasta_files, fasta_id_type, uppercase_symbols) %>%
    # remove decoys from fasta file then drop column
    filter(isdecoy == FALSE) %>%
    select(-isdecoy)

  #### protein-to-accession mapping table
  # convert protein_id to a list of respective accessions  -->>  long-format tibble
  l = strsplit(protein_id, protein_separation_character, fixed = T)
  prot2acc = tibble(
    protein_id = rep(protein_id, lengths(l)),
    acc_input = unlist(l, use.names = F),
    isdecoy = FALSE
  )

  # extract the accession
  prot2acc$acc = fasta_id_short(prot2acc$acc_input, fasta_id_type = fasta_id_type)
  # accession cannot be less than 3 characters
  prot2acc$acc[nchar(prot2acc$acc) < 3] = NA
  # Filter valid accessions by regex, decoy entries are not allowed. eg; REV_ DECOY_
  # Check both full and "short" IDs. eg; suppose that input software prefixes the full ID with reverse: REV_sp|abc123|HUMAN_EXAMPLE -->> the reverse tag is lost in fasta_id_short(...)
  prot2acc$isdecoy = !is.na(prot2acc$acc) & (fasta_identifier_isdecoy(prot2acc$acc_input) | fasta_identifier_isdecoy(prot2acc$acc))

  # finally, remove all invalid mappings and select unique elements (eg; protein-group with same protein accession listed twice is de-duped)
  prot2acc = prot2acc %>%
    filter(!is.na(acc) & isdecoy == FALSE) %>%
    # merge in the fasta metadata
    left_join(fasta, by = c(acc = "idshort")) %>%
    # take unique set, using input with `fasta` sorting as-is = for each idshort, use data from first fasta file (where idshort was found)
    distinct(protein_id, acc_input, .keep_all = T)

  # replace invalid/missing gene symbols with a '-'
  # also removed uniprot-provided gene symbols like MGI:1111. example; ">sp|Q8BR90|CE051_MOUSE UPF0600 protein C5orf51 homolog OS=Mus musculus (Mouse) OX=10090 GN=MGI:2146232 PE=1 SV=1"
  prot2acc$gene[is.na(prot2acc$gene) | nchar(prot2acc$gene) < 2 | grepl("^(MGI|RGD|HGNC):", prot2acc$gene)] = "-"

  # report mapping success rate at accession and protein-group level
  n_acc = nrow(prot2acc)
  n_acc_ok = sum(!is.na(prot2acc$header))
  append_log(sprintf("%d/%d protein accessions and %d/%d protein groups were mapped to provided fasta file(s)",
                     n_acc_ok, n_acc,
                     sum(protein_id %in% (prot2acc %>% filter(!is.na(header)) %>% pull(protein_id))), length(protein_id)), type= "info")
  if(n_acc_ok < 0.5 * n_acc) {
    append_log("More than 50% of all protein accessions could not be mapped to provided fasta file(s). Double-check whether these are the exact same fasta files as used in input dataset's analysis. If fasta files are unknown, use import_fasta() without any parameters to analyze your data without requiring any fasta files", type = "error")
  }
  if(n_acc_ok < 0.9 * n_acc) {
    append_log("More than 10% of protein accessions could not be mapped to provided fasta file(s). Double-check whether these are the exact same fasta files as used in input dataset's analysis", type = "warning")
  }
  if(n_acc_ok < n_acc) {
    print(prot2acc %>% filter(is.na(header)), n=25)
  }


  #### collapse at protein-group level
  tib_result = prot2acc %>%
    group_by(protein_id) %>%
    summarise(
      accessions = paste(acc, collapse = ";"),
      fasta_headers = paste(header, collapse = ";"),
      # full list of gene symbols, '-' where missing, matching number-of-elements in accessions/fasta_headers
      gene_symbols = paste(gene, collapse = ";"),
      gene_symbol_ucount = n_distinct(gene[gene != "-"]),
      # gene_symbols_or_id; take unique genes (remove missing/empty). If there are none, take (all) long-format accessions
      gene_symbols_or_id = ifelse(
        gene_symbol_ucount > 0,
        paste(unique(gene[gene != "-"]), collapse = ";"), # protein_id without gene symbol are ignored / assumed same as gene symbols that we did find
        paste(idlong, collapse = ";")
      )
    )

  rows_fail = is.na(tib_result$gene_symbols_or_id) | tib_result$gene_symbols_or_id == ""
  if(any(rows_fail)) {
    tib_result$fasta_headers[rows_fail] = NA
    tib_result$gene_symbols_or_id[rows_fail] = NA
  }

  # add protein IDs that are in input ('protein_id' parameter of this function) but not in the fasta files and/or excluded by the reverse/contaminant filter above
  pid_miss_input = setdiff(protein_id, tib_result$protein_id)
  if(length(pid_miss_input) > 0) {
    tib_result = bind_rows(tib_result, tibble(protein_id = pid_miss_input, fasta_headers = NA, gene_symbols_or_id = NA))
  }

  rows_miss_fasta = is.na(tib_result$fasta_headers)
  if(any(rows_miss_fasta)) {
    tib_result$fasta_headers[rows_miss_fasta] = tib_result$gene_symbols_or_id[rows_miss_fasta] = tib_result$protein_id[rows_miss_fasta]
  }

  # debug: tib_result %>% filter(gene_symbol_ucount ==0) %>% select(fasta_headers) %>% print(n=100)
  return(tib_result)
}



#' Merge protein identifiers that map to the same unique set of gene symbols
#'
#' @description
#' Suppose there are 3 proteingroups that have the respective unique genes:
#'
#' 1) GRIA2
#' 2) GRIA2
#' 3) GRIA2;GRIA1
#'
#' Then the latter proteingroup is an ambiguous proteingroup where respective peptides were matched to 2 genes
#' (by the upstream software that generated the dataset that you imported with e.g. `import_dataset_diann()`).
#' This function will merge all peptides that belong to proteingroups 1 and 2 into
#' one new proteingroup (that uniquely maps to GRIA2).
#' Proteingroup 3 will remain a distinct proteingroup; only proteingroups with the exact same
#' gene symbol set are matched (technical detail: this uses the definitions in `dataset$proteins$gene_symbols_or_id`).
#'
#' This function will fully update the dataset's protein table and the protein_id column in the peptide table.
#'
#' In most use-cases, this function should be used immediately after `import_fasta()` as shown in the example.
#'
#' @examples \dontrun{
#' # example for using this function. You'll need to update the file paths accordingly
#' library(msdap)
#' dataset = import_dataset_diann("C:/data/diann_report.tsv")
#' dataset = import_fasta(dataset, files = "C:/data/uniprot_human_2020-01.fasta")
#' dataset = merge_proteingroups_by_gene(dataset)
#' }
#' @param dataset a valid dataset. Prior to calling this function, you must import protein metadata from FASTA using the `import_fasta()` function
#' @export
merge_proteingroups_by_gene = function(dataset) {

  # double-check that a fasta has been loaded, and that the dataset object is compatible with msdap 1.2 or later (i.e. has protein data in current format)
  if(!is.data.frame(dataset$proteins) ||
     !all(c("protein_id", "fasta_headers", "gene_symbols", "gene_symbol_ucount", "gene_symbols_or_id") %in% colnames(dataset$proteins)) ||
     !any(is.numeric(dataset$proteins$gene_symbol_ucount) & is.finite(dataset$proteins$gene_symbol_ucount) & dataset$proteins$gene_symbol_ucount > 0)
  ) {
    append_log("dataset is missing essential information in the protein table. No FASTA has been imported for this dataset or the dataset was analyzed with an outdated MS-DAP version.\nTo use this function you won't have to re-run the entire analysis pipeline for this dataset, just run the import_fasta() function to update the current dataset object (assuming this is a dataset searched against a uniprot FASTA)", type = "error")
  }

  if(anyNA(dataset$proteins)) {
    append_log("NA values are not allowed in dataset$proteins", type = "error")
  }


  ### 1) create mapping table from old proteingroup identifiers to new IDs for the updated groups

  protein_to_gene = dataset$peptides %>%
    # count unique peptides per protein_id in dataset$peptides, this dictates the prio/order of uniprot symbols
    distinct(protein_id, peptide_id) %>%
    count(protein_id) %>%
    # after merging proteingroups, the uniprot IDs from the proteingroup with most evidence are in front
    arrange(desc(n)) %>%
    # add gene symbols (using uniprot ID where gene symbols are absent)
    left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by = "protein_id")

  gene_to_new_proteinid__long = protein_to_gene %>%
    select(-n) %>%
    # unlist the uniprot accessions and find set of unique IDs per gene symbol
    # (there may be duplicates across proteingroups, can't trust input data)
    mutate(protein_id_new = strsplit(protein_id, " *; *")) %>%
    tidyr::unchop(protein_id_new) %>%
    distinct(gene_symbols_or_id, protein_id_new)

  gene_to_new_proteinid = gene_to_new_proteinid__long %>%
    # collapse into a string: these are the new protein_id for respective gene_symbols_or_id
    group_by(gene_symbols_or_id) %>%
    summarise(protein_id_new = paste(protein_id_new, collapse = ";"), .groups = "drop")

  # finally, merge both tables and print stats to console
  protein_to_new_proteinid = left_join(protein_to_gene, gene_to_new_proteinid, by = "gene_symbols_or_id")

  tmp = sum(protein_to_new_proteinid$protein_id != protein_to_new_proteinid$protein_id_new)
  append_log(sprintf("%d / %d (%.1f%%) proteingroups that matched to the same gene symbol(s) were merged",
                     tmp, nrow(protein_to_new_proteinid), tmp / nrow(protein_to_new_proteinid) * 100), type = "info")


  ### 2) update dataset$peptides : replace protein_id

  i = data.table::chmatch(dataset$peptides$protein_id, protein_to_new_proteinid$protein_id)
  stopifnot(!is.na(i))
  dataset$peptides$protein_id = protein_to_new_proteinid$protein_id_new[i]


  ### 3) update dataset$proteins : deconstruct and merge protein data

  ## example data for a proteingroup with multiple accessions, and some of which have no gene symbol
  # $ protein_id        : chr "A0A087WUQ1;B3KPU0;E5RGS7;H7C1D1"
  # $ accessions        : chr "A0A087WUQ1;B3KPU0;E5RGS7;H7C1D1"
  # $ fasta_headers     : chr ">tr|A0A087WUQ1|A0A087WUQ1_HUMAN Isoform of Q96JY6, PDZ and LIM domain 2 OS=Homo sapiens OX=9606 GN=PDLIM2 PE=1 "| __truncated__
  # $ gene_symbols      : chr "PDLIM2;PDLIM2;-;-"
  # $ gene_symbol_ucount: int 1
  # $ gene_symbols_or_id: chr "PDLIM2"

  # first, unpack the current protein data. These columns are semicolon-delimited and assumed to have same length
  p = dataset$proteins %>%
    select(protein_id, fasta_headers, gene_symbols) %>%
    mutate(
      protein_id = strsplit(protein_id, ";", fixed = TRUE),
      fasta_headers = strsplit(fasta_headers, ";", fixed = TRUE),
      gene_symbols = strsplit(gene_symbols, ";", fixed = TRUE)
    )
  stopifnot(lengths(p$protein_id) == lengths(p$fasta_headers))
  stopifnot(lengths(p$protein_id) == lengths(p$gene_symbols))

  p = p %>%
    tidyr::unchop(cols = c(protein_id, fasta_headers, gene_symbols)) %>%
    # importantly, take unique: the same uniprot ID may be part of multiple proteingroups
    distinct(protein_id, .keep_all = TRUE)

  # join the fasta/symbol metadata to the new proteingroup identifiers
  newproteins = gene_to_new_proteinid__long %>%
    rename(protein_id = protein_id_new) %>%
    left_join(p, by = "protein_id") %>%
    group_by(gene_symbols_or_id) %>%
    summarise(
      fasta_headers = paste(fasta_headers, collapse = ";"),
      gene_symbols = paste(gene_symbols, collapse = ";"),
      ## definition of ucount analogous to import_protein_metadata_from_fasta()
      gene_symbol_ucount = n_distinct(gene_symbols[gene_symbols != "-"])
    )

  dataset$proteins = gene_to_new_proteinid %>%
    rename(protein_id = protein_id_new) %>%
    mutate(accessions = protein_id) %>%
    left_join(newproteins, by = "gene_symbols_or_id") %>%
    select(protein_id, accessions, fasta_headers, gene_symbols, gene_symbol_ucount, gene_symbols_or_id)

  # QC: no NA's (i.e. no issues while joining tables)  +  all updated protein_id in peptide table are present in protein table
  stopifnot(!anyNA(dataset$proteins))
  stopifnot(unique(dataset$peptides$protein_id) %in% dataset$proteins$protein_id)

  return(dataset)
}



#' merge fractionated samples
#'
#' to collapse sample fractions, the sample metadata table must contain columns 'fraction' and 'shortname' (eg; sample_id is unique sample/measurement, whereas shortname is the identifier to which respective fractions are collapsed to)
#' @param dataset a valid dataset
#' @importFrom data.table chmatch
#' @export
merge_fractionated_samples = function(dataset) {
  if(!all(c("shortname", "fraction") %in% colnames(dataset$samples))) {
    append_log("to collapse sample fractions, the sample metadata table must contain columns 'fraction' and 'shortname' (eg; sample_id is unique sample/measurement, whereas shortname is the identifier to which respective fractions are collapsed to)", type = "error")
  }


  # prior to replacing sample_id, store mappings
  sid_shortname = dataset$samples %>% distinct(sample_id, shortname)

  # 1) replace sample_id by shortname in both the samples and peptides tibbles
  dataset$peptides$sample_id = dataset$samples$shortname[data.table::chmatch(dataset$peptides$sample_id, dataset$samples$sample_id)]
  dataset$samples = dataset$samples %>% mutate(sample_id = shortname) %>% select(-fraction) %>% distinct(sample_id, .keep_all = T)

  # 2) collapse by sample_id*peptide_id
  dataset$peptides = peptides_collapse_by_sequence(dataset$peptides, prop_peptide = "peptide_id")

  # 3) update peptide/protein counts per sample and remove any caching, if present
  dataset$samples = peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))
  dataset = invalidate_cache(dataset)

  # 4) if any contrasts are present, update the respective sample identifiers
  if(is.list(dataset$contrasts) && length(dataset$contrasts) > 0) {
    for(index in seq_along(dataset$contrasts)) {
      contr = dataset$contrasts[[index]]
      # replace sample_id with shortname
      contr$sampleid_condition1 = unique(sid_shortname$shortname[match(contr$sampleid_condition1, sid_shortname$sample_id)])
      contr$sampleid_condition2 = unique(sid_shortname$shortname[match(contr$sampleid_condition2, sid_shortname$sample_id)])
      contr$sample_table = contr$sample_table %>% mutate(sample_id = sid_shortname$shortname[match(sample_id, sid_shortname$sample_id)])
      # keep track of removed samples (duplicates) @ `sample_table` --> sync with `model_matrix` which should be aligned
      row_isdupe = duplicated(contr$sample_table$sample_id)
      contr$sample_table = contr$sample_table[!row_isdupe,]
      contr$model_matrix = contr$model_matrix[!row_isdupe,]
      if(!is_tibble(contr$sample_table)) { # setting rownames on a tibble is deprecated, but we can/should on data.frame
        rownames(contr$sample_table) = seq_len(nrow(contr$sample_table))
      }
      rownames(contr$model_matrix) = seq_len(nrow(contr$sample_table))
      # update contrast
      dataset$contrasts[[index]] = contr
    }
  }


  return(dataset)
}



#' highly optimized aggregation of the peptides data table.
#' rewrite aggregate_tibble_by_datatables; more speed, support charge/mz columns (if present), take RT from most abundant peptide with prio for detected. the latter matters when merging data for plain sequences where one precursor has an MS/MS identification and the other is MBR
#'
#' @param peptides todo
#' @param prop_peptide must be either 'sequence_plain', 'sequence_modified' or 'peptide_id' (the former 2 are the common use-case)
#'
#' @importFrom data.table data.table chmatch setkey setorder
aggregate_tibble_by_datatables = function(peptides, prop_peptide = "sequence_plain") {
  if(!(prop_peptide %in% c("sequence_plain", "sequence_modified", "peptide_id"))) {
    append_log("prop_peptide parameter must be either 'sequence_plain', 'sequence_modified' or 'peptide_id' (the former 2 are the common use-case)", type = "error")
  }

  # convert input tibble/data.frame to data.table and create our own numeric key (for re-use later)
  x = data.table::data.table(peptides, stringsAsFactors = F)
  data.table::setkeyv(x, c("sample_id", prop_peptide)) # data.table::key(x)
  x[ , DTkeyindex := .GRP, by=eval(c("sample_id", prop_peptide))]
  ## previous implementation, not using data.table built-in functions, slightly slower
  # x$DTkeystring = paste0(x$sample_id,"#",x[[prop_peptide]])
  # x$DTkeyindex2 = data.table::chmatch(x$DTkeystring, unique(x$DTkeystring))
  # data.table::setkey(x, "DTkeyindex2")
  # all(x$DTkeyindex == x$DTkeyindex2) # check that this is the same as new group keys above

  # transform intensities, then collapse each peptide*sample; check if any detect + sum intensities
  x$intensity = 2^(x$intensity)
  y = x[ , .(detect = any(detect, na.rm=T), intensity = sum(intensity, na.rm=T)), by = DTkeyindex]

  ## sort table and add first encounter of some value
  # detect prio, then intensity descending -> first RT.  take most abundant peak for reference RT of peptide p in sample s instead of their mean, just in case these are far apart (in which case the mean RT makes no sense)
  data.table::setorder(x, -detect, -intensity, na.last = T) # example code to double-check sorting: x = data.table(a=1:3, b=c(F,T,F));x;setorder(x, -b,-a, na.last=FALSE);x
  i = match(y$DTkeyindex, x$DTkeyindex)
  y$rt = x$rt[i]
  y$protein_id = x$protein_id[i]
  y$sample_id = x$sample_id[i]

  if (prop_peptide == "peptide_id" && "charge" %in% colnames(x)) {
    y$charge = x$charge[i]
  }

  if("mz" %in% colnames(x)) {
    y$mz = x$mz[i]
  }

  if("isdecoy" %in% colnames(x)) {
    y$isdecoy = x$isdecoy[i]
  } else {
    y$isdecoy = FALSE
  }

  # set the peptide_id to peptide identity used for grouping
  y$peptide_id = x[[prop_peptide]][i]
  # plain sequence always same as input
  y$sequence_plain = x$sequence_plain[i]
  # if grouped by plain sequence, there is no modified sequence anymore
  if (prop_peptide == "sequence_plain") {
    y$sequence_modified = y$sequence_plain
  } else {
    y$sequence_modified = x$sequence_modified[i]
  }

  ## !!! careful, potentially re-ordering the data.table below so index lookup i is no longer valid from here on
  rm(i)
  if("confidence" %in% colnames(x)) {
    # confidence ascending -> take first value
    data.table::setorder(x, confidence, na.last = T)
    y$confidence = x$confidence[match(y$DTkeyindex, x$DTkeyindex)]
  } else {
    y$confidence = NA
  }

  # analogous
  if("confidence_score" %in% colnames(x)) {
    # confidence_score descending -> take first value
    data.table::setorder(x, -confidence_score, na.last = T) # note the - prefix for descending
    y$confidence_score = x$confidence_score[match(y$DTkeyindex, x$DTkeyindex)]
  } # not a "standard field" so don't have to append NA column if missing


  # stopifnot(anyDuplicated(y$DTkeyindex) == 0)
  # stopifnot(anyDuplicated(y[,c("sequence_modified","sample_id")]) == 0)

  ## enforce data integrity (detect must be T/F, RT has to be finite or NA)
  # y$rt[!is.finite(y$rt)] = NA
  y$detect = y$detect %in% TRUE # replace NA with false. faster than checking !is.na(x) & x==T. example; c(T,F,0,NA,1) %in% TRUE
  # tranform intensity back to log2
  y$intensity = log2(y$intensity)
  # drop sample*peptide lookup key from result table
  y$DTkeyindex = NULL
  return(y)
}


# select_topn_peptide_per_protein_legacy = function(peptides, samples, topn_peptides = 10, topn_all_groups = F) {
#   if(!is.finite(topn_peptides) || topn_peptides < 1) {
#     return(peptides)
#   }
#
#   tib = peptides %>% left_join(samples %>% select(sample_id, group) %>% add_count(group, name = "group_size"), by = "sample_id")
#   # filter rules
#   ngrp = length(unique(tib$group))
#   ngrp_topn = ifelse(topn_all_groups, ngrp, ceiling(ngrp * 0.5))
#
#   #### score: peptide*group
#   # score components; abundance (median intensity) and reproducibility (number of detects  &  CoV)
#   tib_peptide_group = tib %>%
#     group_by(peptide_id, group) %>%
#     summarise(
#       score_intensity = stats::median(intensity),
#       score_n = sum(detect) / group_size[1], # fraction of samples where detected
#       #score_n = n() / group_size[1],
#       score_cov = coefficient_of_variation(2^intensity)
#     )
#
#   # overall score distributions
#   distr_intensity = ecdf(peptides$intensity) # overall peptide intensities
#   distr_cov = ecdf(log10(tib_peptide_group$score_cov)) # CoV's from each peptide * group
#   # plot(distr_intensity); plot(distr_cov)
#
#   # scale scores using distributions, where smaller cov value = better score
#   tib_peptide_group$score_intensity = distr_intensity(tib_peptide_group$score_intensity)
#   tib_peptide_group$score_cov = 1 - distr_cov(log10(tib_peptide_group$score_cov))
#
#   # overall score = average of 3 score components
#   tib_peptide_group$score_n[!is.finite(tib_peptide_group$score_n)] = 0
#   tib_peptide_group$score_cov[!is.finite(tib_peptide_group$score_cov)] = 0
#   tib_peptide_group$score_intensity[!is.finite(tib_peptide_group$score_intensity)] = 0
#   # simply use average over 'fraction of replicated where quantified', CoV and intensity
#   tib_peptide_group$score = (tib_peptide_group$score_n + tib_peptide_group$score_cov + tib_peptide_group$score_intensity) / 3
#   # hist(tib_peptide_group$score)
#
#
#   #### score: peptide
#   # select topn groups for each peptide (given score attribute)
#   tib_peptide_group_topn = tib_peptide_group %>%
#     group_by(peptide_id) %>%
#     top_n(ngrp_topn, wt = score)
#   # if multiple groups are to be used for computing peptide score, summarize by taking the mean score
#   if (ngrp_topn > 1) {
#     tib_peptide_group_topn = tib_peptide_group_topn %>% summarise(score = mean(score, na.rm = T)) # already grouped by peptide_id
#   }
#
#   #### score: foreach protein, select topN peptides
#   tib_peptide_group_topn$protein_id = peptides$protein_id[match(tib_peptide_group_topn$peptide_id, peptides$peptide_id)]
#   tib_valid_peptides = tib_peptide_group_topn %>%
#     group_by(protein_id) %>%
#     top_n(topn_peptides, wt = score)
#
#   return(tib_valid_peptides$peptide_id)
#   # #### finally, return subset of input peptides that remain in the 'topn filtered peptide tibble'
#   # peptides %>% filter(peptide_id %in% tib_valid_peptides$peptide_id)
# }
