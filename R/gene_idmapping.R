
#' Parse HGNC gene identifier lookup table that was downloaded from genenames.org
#'
#' download link: https://www.genenames.org/download/statistics-and-files/
#' table: "Complete dataset download links" -->> "Complete HGNC approved dataset text json" -->> download the "TXT" table
#' filename is typically something like hgnc_complete_set.txt
#' URL as of September 2023; https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
#'
#' # alternatively;
#' table: "Total Approved Symbols" -->> "TXT" / "text file in TSV format"
#' filename is typically something like non_alt_loci_set.txt
#'
#' @param f full path to the downloaded table (expected to be tsv format)
#' @returns a long-format table with columns; hgnc_id, hgnc_symbol, type, value
#' @export
hgnc_lookuptable = function(f) {
  check_parameter_is_string(f)
  if(!file.exists(f)) {
    append_log(paste("file does not exist:", f), type = "error")
  }

  result = NULL

  # parse HGNC table from disk
  hgnc = data.table::fread(f, data.table = F, stringsAsFactors = F)

  # check that all required columns are present + rename columns
  cols = list(
    hgnc_id = c("hgnc_id", "HGNC ID"),
    hgnc_symbol = c("symbol", "Approved symbol"),
    alias_symbol = c("alias_symbol", "Previous symbols"),
    prev_symbol = c("prev_symbol", "Alias symbols"),
    ensembl_id = c("ensembl_id", "Ensembl gene ID", "Ensembl ID(supplied by Ensembl)"),
    entrez_id = c("entrez_id", "NCBI Gene ID", "NCBI Gene ID(supplied by NCBI)"),
    mgi_id = c("mgi_id", "mgd_id", "Mouse genome database ID", "Mouse genome database ID(supplied by MGI)"),
    rgd_id = c("rgd_id", "Rat genome database ID", "Rat genome database ID(supplied by RGD)")
  )
  hgnc = robust_header_matching(hgnc, column_spec = cols, columns_required = c("hgnc_id", "hgnc_symbol"))

  hgnc$hgnc_symbol = toupper(hgnc$hgnc_symbol)

  # remove invalid rows (this shouldn't occur but check anyway)
  hgnc = hgnc %>% filter(!is.na(hgnc_id) & hgnc_id != "" & !is.na(hgnc_symbol) & nchar(hgnc_symbol) >= 2)

  # no optional columns, done
  if(ncol(hgnc) == 2) {
    return(hgnc)
  }

  if(any(c("alias_symbol", "prev_symbol") %in% colnames(hgnc))) {
    # collect synonyms from available columns
    if("alias_symbol" %in% colnames(hgnc)) {
      hgnc$synonyms = replace(hgnc$alias_symbol, list = is.na(hgnc$alias_symbol), values = "")
      if("prev_symbol" %in% colnames(hgnc)) {
        hgnc$synonyms = paste(replace(hgnc$prev_symbol, list = is.na(hgnc$prev_symbol), values = ""), hgnc$synonyms)
      }
    } else {
      hgnc$synonyms = hgnc$prev_symbol
    }

    # parse synonyms into a long-format lookup table & remove empty / ambiguous
    # (synonyms that overlap with main symbol OR are found 2+ times)
    l = strsplit(toupper(hgnc$synonyms), "[ ,;|]+") # support various delimiters
    map_synonym = tibble::tibble(
      hgnc_id = rep(hgnc$hgnc_id, lengths(l)),
      symbol = unlist(l, recursive = F, use.names = F)
    ) %>%
      dplyr::filter(nchar(symbol) >= 2 &  ! symbol %in% hgnc$hgnc_symbol) %>%
      # deal with duplicate synonyms within 1 row/hgnc_id
      # (e.g. due to our case conversion or overlap between previous and alias columns)
      dplyr::distinct_all()

    # count how often each synonym is found --> remove those with multiple entries
    symbol_ambiguous = map_synonym %>% dplyr::count(symbol) %>% dplyr::filter(n > 1) %>% pull(symbol)
    map_synonym = map_synonym %>% dplyr::filter( ! symbol %in% symbol_ambiguous)

    result = bind_rows(
      result,
      map_synonym %>% mutate(type = "synonym") %>% select(hgnc_id, type, value = symbol)
    )
  }



  # format/filter crossreference ids
  if("ensembl_id" %in% colnames(hgnc)) {
    hgnc$ensembl_id = gsub("[ ,;/|].*", "", hgnc$ensembl_id) # remove multiple mappings
    hgnc$ensembl_id[is.na(hgnc$ensembl_id) | nchar(hgnc$ensembl_id) < 5] = NA # invalid IDs to NA
    result = bind_rows(
      result,
      hgnc %>% filter(!is.na(ensembl_id)) %>% mutate(type = "ensembl_id") %>% select(hgnc_id, type, value = ensembl_id)
    )
  }
  if("entrez_id" %in% colnames(hgnc)) {
    hgnc$entrez_id = gsub("\\D.*", "", as.character(hgnc$entrez_id)) # remove multiple mappings
    hgnc$entrez_id[is.na(hgnc$entrez_id) | !grepl("^\\d+$", hgnc$entrez_id)] = NA # invalid IDs to NA
    result = bind_rows(
      result,
      hgnc %>% filter(!is.na(entrez_id)) %>% mutate(type = "entrez_id") %>% select(hgnc_id, type, value = entrez_id)
    )
  }
  # strip and reconstruct 'MGI:' prefix to enforce + convert 1:n mappings to long format via string-split
  if("mgi_id" %in% colnames(hgnc)) {
    l = strsplit(gsub("MG[ID]: *", "", hgnc$mgi_id, ignore.case = T), "[ ,;/|]+")
    tmp = tibble::tibble(
      hgnc_id = rep(hgnc$hgnc_id, lengths(l)),
      id = unlist(l, recursive = F, use.names = F)
    ) %>%
      dplyr::filter(!is.na(id) & id != "") %>%
      dplyr::distinct_all() %>%
      mutate(id = paste0("MGI:", id))

    result = bind_rows(
      result,
      tmp %>% mutate(type = "mgi_id") %>% select(hgnc_id, type, value = id)
    )
  }
  # analogous to mgi_id
  if("rgd_id" %in% colnames(hgnc)) {
    l = strsplit(gsub("RGD: *", "", hgnc$rgd_id, ignore.case = T), "[ ,;/|]+")
    tmp = tibble::tibble(
      hgnc_id = rep(hgnc$hgnc_id, lengths(l)),
      id = unlist(l, recursive = F, use.names = F)
    ) %>%
      dplyr::filter(!is.na(id) & id != "") %>%
      dplyr::distinct_all() %>%
      mutate(id = paste0("RGD:", id))

    result = bind_rows(
      result,
      tmp %>% mutate(type = "rgd_id") %>% select(hgnc_id, type, value = id)
    )
  }

  if(!is.data.frame(result) || nrow(result) == 0) {
    return(hgnc %>% select(hgnc_id, hgnc_symbol))
  }

  # add hgnc_id that do not have any synonyms or mappings to mgi/rgd/entrez/ensembl
  id_missing = setdiff(hgnc$hgnc_id, result$hgnc_id)
  if(length(id_missing) > 0) {
    id_missing__symbol = hgnc$hgnc_symbol[match(id_missing, hgnc$hgnc_id)]
    result = bind_rows(
      result,
      tibble::tibble(hgnc_id = id_missing, type = "onlysymbol", value = id_missing__symbol)
    )
  }

  # add HGNC symbols
  result$hgnc_symbol = hgnc$hgnc_symbol[match(result$hgnc_id, hgnc$hgnc_id)]
  # return results with re-ordered column names
  result %>% select(hgnc_id, hgnc_symbol, type, value)
}



#' Parse MGI ID mapping table table into a lookup table with symbols and uniprot IDs
#'
#' download link: https://www.informatics.jax.org/downloads/reports/index.html
#' table: "MGI Marker associations to SWISS-PROT and TrEMBL protein IDs (tab-delimited)"
#' filename is typically something like MRK_SwissProt_TrEMBL.rpt
#'
#' @param f full path to the downloaded table (expected to be tsv format)
#' @export
mgi_lookuptable = function(f) {
  check_parameter_is_string(f)
  if(!file.exists(f)) {
    append_log(paste("file does not exist:", f), type = "error")
  }

  # MGI Marker Accession ID | Marker Symbol | Status | Marker Name | cM Position | Chromosome | SWISS-PROT/TrEMBL Protein Accession IDs (space-delimited)
  x = data.table::fread(f, data.table = FALSE, stringsAsFactors = FALSE)
  if(ncol(x) != 7) {
    append_log(paste0("invalid MGI table; we expected 7 columns but found ", ncol(x), ". Is this MRK_SwissProt_TrEMBL.rpt @ MGI Marker associations to SWISS-PROT and TrEMBL protein IDs (tab-delimited)  ?"), type = "error")
  }

  colnames(x) = c("mgi_id", "symbol", "status", "marker_name", "cm_position", "chromosome", "uniprot_id")

  # enforce MGI: prefix for identifiers
  x$mgi_id = paste0("MGI:", gsub("MGI: *", "", x$mgi_id, ignore.case = T))
  x$symbol = toupper(x$symbol)
  x$uniprot_id = strsplit(x$uniprot_id, "[ ,;]+")
  y = x %>%
    select(mgi_id, symbol, uniprot_id) %>%
    tibble::as_tibble() %>%
    tidyr::unchop(uniprot_id, keep_empty = TRUE)
  y$uniprot_id[is.na(y$uniprot_id)] = ""
  return(y)
}



#' Parse RGD ID mapping table table into a lookup table with symbols and uniprot IDs
#'
#' download link: https://download.rgd.mcw.edu/data_release/GENES_RAT.txt
#' README: https://download.rgd.mcw.edu/data_release/GENES_README
#'
#' @param f full path to the downloaded table (expected to be tsv format)
#' @export
rgd_lookuptable = function(f) {
  check_parameter_is_string(f)
  if(!file.exists(f)) {
    append_log(paste("file does not exist:", f), type = "error")
  }

  # load table from disk
  x = data.table::fread(f, data.table = FALSE, stringsAsFactors = FALSE)

  # check that all required columns are present + rename columns
  cols = list(
    rgd_id = c("GENE_RGD_ID"),
    symbol = c("symbol", "SYMBOL"),
    uniprot_id = c("UNIPROT_ID")
  )
  x = robust_header_matching(x, column_spec = cols, columns_required = names(cols))

  # enforce RGD: prefix for identifiers
  x$rgd_id = paste0("RGD:", gsub("\\D+", "", x$rgd_id, ignore.case = T))
  x$symbol = toupper(x$symbol)
  x$uniprot_id = strsplit(x$uniprot_id, "[ ,;]+")
  y = x %>%
    tibble::as_tibble() %>%
    tidyr::unchop(uniprot_id, keep_empty = TRUE)
  y$uniprot_id[is.na(y$uniprot_id)] = ""
  return(y)
}



#' Map the protein_id column, which are assumed to contain uniprot mouse/rat identifiers, in a table to HGNC human gene IDs
#'
#' This function is currently limited such that the first matching human gene per proteingroup is returned,
#' so ambiguous proteingroups with multiple gene symbols (e.g. "GRIA1;GRIA2") might not always yield the 'leading'
#' gene and additional/ambiguous genes are disregarded). Columns describing the HGNC (genenames.org) ID,
#' NCBI Entez ID and Ensembl gene ID are appended to the provided table.
#'
#' Workflow for ID mapping, where subsequent steps are fallbacks;
#' 1) uniprot id -> mgi/rgd id -> hgnc id
#' 2) uniprot symbol -> mgi/rgd id -> hgnc id
#' 3) uniprot symbol -> hgnc official symbol -> hgnc id
#' 4) uniprot symbol -> hgnc non-ambiguous synonyms -> hgnc id
#'
#' @examples \dontrun{
#'   # prior to this, import a dataset/samples/fasta & setup contrasts & apply filter_datasets()
#'   # alternatively, use analysis_quickstart() to do everything
#'   dataset = dea(dataset, dea_algorithm = "deqms")
#'   dataset = differential_detect(dataset, min_peptides_observed = 2, min_samples_observed = 3,
#'                                 min_fraction_observed = 0.5, return_wide_format = FALSE)
#'   print_dataset_summary(dataset)
#'
#'   x = summarise_stats(dataset, return_dea = TRUE, return_diffdetect = TRUE,
#'                       diffdetect_zscore_threshold = 5, remove_ambiguous_proteingroups = TRUE)
#'   hgnc = hgnc_lookuptable(f = "C:/DATA/hgnc/hgnc__non_alt_loci_set___2023-03-21.tsv")
#'   mgi = mgi_lookuptable(f = "C:/DATA/mgi/2023-07-15/MRK_SwissProt_TrEMBL.rpt")
#'   y = protein2gene_orthologs(x, hgnc, mgi)
#' }
#' @param x a data.table with columns protein_id and symbol (e.g. the "leading symbol" extracted from the uniprot gene symbols @ dataset$proteins$gene_symbols_or_id)
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
#' @param idmap MGI lookup table from `mgi_lookuptable()` or a RGD lookup table from `rgd_lookuptable()`
#' @param remove_nohgnc set to TRUE to remove all rows that could not be matched to a HGNC ID. default; FALSE
#' @export
protein2gene_orthologs = function(x, hgnc, idmap, remove_nohgnc = FALSE) {
  if(!is.data.frame(x) || nrow(x) == 0 || !all(c("protein_id", "symbol") %in% colnames(x)) ) {
    append_log("x must be a data.frame with columns 'protein_id' and 'symbol'", type = "error")
  }
  if(anyNA(x$protein_id) || !is.character(x$protein_id) || any(x$protein_id == "")) {
    append_log("column 'protein_id' in x must be of character type and not contain any NA values or empty strings", type = "error")
  }
  if(!is.data.frame(hgnc) || nrow(hgnc) == 0 || !all(c("hgnc_id", "hgnc_symbol", "type", "value") %in% colnames(hgnc)) ) {
    append_log("hgnc must be a data.frame with columns 'hgnc_id', 'hgnc_symbol', 'type', 'value' as typically prepared using the hgnc_lookuptable() function", type = "error")
  }
  if(!is.data.frame(idmap) || nrow(idmap) == 0 || !all(c("symbol", "uniprot_id") %in% colnames(idmap)) || !colnames(idmap)[1] %in% hgnc$type) {
    append_log("idmap must be a data.frame with columns 'symbol', 'uniprot_id' and the first column should contain identifiers that are also in the hgnc table (e.g. first column name should be 'mgi_id' or 'rgd_id'). e.g. as obtained from the mgi_lookuptable() function", type = "error")
  }

  # remove preexisting columns, if any. e.g. suppose that `x` contains entrez_id column but here we don't have entrez_id
  # in the hgnc table; this table never gets overwritten
  x$hgnc_id = x$hgnc_symbol = x$ensembl_id = x$entrez_id = NULL

  # replace empty symbols, if any
  x$symbol_input = x$symbol
  rows_symbol_fail = is.na(x$symbol) | nchar(x$symbol) < 2
  x$symbol[rows_symbol_fail] = x$protein_id[rows_symbol_fail]
  x$symbol = toupper(x$symbol) # enforce uppercase

  proteins_long = x %>%
    select(protein_id, symbol) %>%
    distinct(protein_id, .keep_all = TRUE) %>%
    mutate(uniprot_id = strsplit(protein_id, split = " *; *")) %>%
    tidyr::unchop(uniprot_id, keep_empty = TRUE) %>%
    filter(nchar(uniprot_id) >= 3) %>% # shouldn't be needed but enforce check for empty IDs anyway
    mutate(
      uniprot_id = gsub("-.*", "", uniprot_id) #,
      # by_symbol = TRUE # default = leftover matching by symbol
    )

  ## map from uniprot ID to MGI/RGD
  i = match(proteins_long$uniprot_id, idmap$uniprot_id)              # match by uniprot id
  i[is.na(i)] = match(proteins_long$symbol[is.na(i)], toupper(idmap$symbol))  # fallback matching by symbol (should align well between MGI and uniprot mouse proteome)
  proteins_long$xref_id = unlist(idmap[i,1], recursive = F, use.names = F) # first column in the data.frame / tibble

  ## map to HGNC
  # hgnc subset of relevant crossreference identifiers
  hgnc_xref = hgnc %>% filter(type == colnames(idmap)[1])
  proteins_long$hgnc_id = hgnc_xref$hgnc_id[match(proteins_long$xref_id, hgnc_xref$value)]
  # flag all protein_id that we directly matched as 'not matched by symbols but rather through identifiers'
  # proteins_long$by_symbol[!is.na(proteins_long$hgnc_id)] = FALSE

  # fallback matching by official gene symbol
  # note that we shouldn't use subsets of the hgnc table, e.g. for crossreference or synonyms, as not all HGNC entries might have crossreferences/synonyms
  proteins_long$hgnc_id[is.na(proteins_long$hgnc_id)] = hgnc$hgnc_id[match(proteins_long$symbol[is.na(proteins_long$hgnc_id)], hgnc$hgnc_symbol)]

  # fallback matching by synonym
  rows = is.na(proteins_long$hgnc_id)
  if(any(rows)) {
    hgnc_synonym = hgnc %>% filter(type == "synonym")
    if(nrow(hgnc_synonym) > 0) {
      proteins_long$hgnc_id[rows] = hgnc_synonym$hgnc_id[match(proteins_long$symbol[rows], hgnc_synonym$value)]
    }
  }


  ### mapping table from proteingroups (protein_id) to HGNC
  # importantly, we don't deal with ambiguous mappings atm
  # e.g. if a proteingroup contains multiple gene symbols, we just return the first human gene ID we found
  protein_to_gene = proteins_long %>%
    filter(!is.na(hgnc_id)) %>%
    # retain just 1 entry per proteingroup -->> any 1:n mapping is lost
    distinct(protein_id, .keep_all = TRUE)
  # direct 1:1 matching
  i = match(x$protein_id, protein_to_gene$protein_id)
  x$hgnc_id = protein_to_gene$hgnc_id[i]

  ### add hgnc_symbol and entrez/ensembl IDs if present in hgnc table
  x = hgnc_add_xrefs(x, hgnc)

  ### restore input symbols
  x$symbol = x$symbol_input
  x$symbol_input = NULL

  # log status to console. use distinct(protein_id) because some input tables might contain multiple entries for
  # some protein_id, e.g. dataset$de_proteins is provided as-is when it contains 2 contrasts
  tmp = x %>% distinct(protein_id, .keep_all = TRUE)
  n_input = nrow(tmp)
  n_fail = sum(is.na(tmp$hgnc_id))
  cat(sprintf("%d / %d (%.1f%%) unique protein_id could not be mapped to a HGNC human gene ID\n",
              n_fail, n_input, n_fail / n_input * 100 ))

  if(remove_nohgnc) {
    x = x %>% filter(!is.na(hgnc_id))
  }
  return(x)
}



#' Map the the symbol column in a table to HGNC human gene IDs by matching official gene symbols and synonyms
#'
#' @param x a data.table with a column symbol
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
#' @export
protein2gene_by_symbol = function(x, hgnc) {
  if(!is.data.frame(x) || nrow(x) == 0 || !all("symbol" %in% colnames(x) || !is.character(x$symbol)) ) {
    append_log("x must be a data.frame with column 'symbol' (character type)", type = "error")
  }
  if(!is.data.frame(hgnc) || nrow(hgnc) == 0 || !all(c("hgnc_id", "hgnc_symbol", "type", "value") %in% colnames(hgnc)) ) {
    append_log("hgnc must be a data.frame with columns 'hgnc_id', 'hgnc_symbol', 'type', 'value' as typically prepared using the hgnc_lookuptable() function", type = "error")
  }

  x$symbol_input = x$symbol
  rows_symbol_fail = is.na(x$symbol) | nchar(x$symbol) < 2
  x$symbol = toupper(x$symbol)

  # match directly by official gene symbol
  i = match(x$symbol, hgnc$hgnc_symbol)
  x$hgnc_id = hgnc$hgnc_id[i]

  # for rows that failed to match, try matching against the HGNC synonym table
  rows = is.na(x$hgnc_id)
  if(any(rows)) {
    hgnc_synonym = hgnc %>% filter(type == "synonym")
    if(nrow(hgnc_synonym) > 0) {
      i = match(x$symbol[rows], hgnc_synonym$value)
      x$hgnc_id[rows] = hgnc_synonym$hgnc_id[i]
    }
  }

  # erase whatever happened to invalid input
  x$hgnc_id[rows_symbol_fail] = NA

  # add hgnc_symbol and entrez/ensembl IDs if present in hgnc table
  x = hgnc_add_xrefs(x, hgnc)

  # restore input symbols
  x$symbol = x$symbol_input
  x$symbol_input = NULL

  # log status to console. use distinct(protein_id) because some input tables might contain multiple entries for
  # some protein_id, e.g. dataset$de_proteins is provided as-is when it contains 2 contrasts
  tmp = x %>% filter(rows_symbol_fail == FALSE) %>% distinct(symbol, .keep_all = TRUE)
  n_input = nrow(tmp)
  n_fail = sum(is.na(tmp$hgnc_id))
  cat(sprintf("%d / %d (%.1f%%) unique symbols could not be mapped to a HGNC human gene ID\n",
              n_fail, n_input, n_fail / n_input * 100 ))

  return(x)
}



#' add hgnc_symbol and entrez/ensembl IDs if present in hgnc table
#'
#' @param x a data.table with column hgnc_id
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
hgnc_add_xrefs = function(x, hgnc) {
  # symbol
  i = match(x$hgnc_id, hgnc$hgnc_id)
  x$hgnc_symbol = hgnc$hgnc_symbol[i]

  # remove optional return values if present
  x$entrez_id = NULL
  x$ensembl_id = NULL

  # subset HGNC lookup table by respective crossref IDs
  hgnc_entrez = hgnc %>% filter(type == "entrez_id")
  hgnc_ensembl = hgnc %>% filter(type == "ensembl_id")

  if(nrow(hgnc_entrez) > 0) {
    i = match(x$hgnc_id, hgnc_entrez$hgnc_id)
    x$entrez_id = hgnc_entrez$value[i]
  }
  if(nrow(hgnc_ensembl) > 0) {
    i = match(x$hgnc_id, hgnc_ensembl$hgnc_id)
    x$ensembl_id = hgnc_ensembl$value[i]
  }

  return(x)
}
