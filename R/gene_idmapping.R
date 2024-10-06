
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
#' @param uppercase_symbols convert all gene symbols to upper case? default: TRUE
#' @returns a long-format table with columns; hgnc_id, hgnc_symbol, type, value
#' @export
hgnc_lookuptable = function(f, uppercase_symbols = TRUE) {
  result = NULL
  check_parameter_is_string(f)
  check_parameter_is_boolean(uppercase_symbols)

  # parse HGNC table from disk
  f = path_exists(f, NULL, try_compressed = TRUE)
  hgnc = as_tibble(read_textfile_compressed(f, as_table = TRUE))
  # if(!file.exists(f)) {
  #   append_log(paste("file does not exist:", f), type = "error")
  # }
  # hgnc = data.table::fread(f, data.table = F, stringsAsFactors = F)


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

  if(uppercase_symbols) {
    hgnc$hgnc_symbol = toupper(hgnc$hgnc_symbol)
  }

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

  # parse HGNC table from disk
  f = path_exists(f, NULL, try_compressed = TRUE)
  x = as_tibble(read_textfile_compressed(f, as_table = TRUE, header = FALSE))
  # if(!file.exists(f)) {
  #   append_log(paste("file does not exist:", f), type = "error")
  # }
  # x = data.table::fread(f, data.table = FALSE, stringsAsFactors = FALSE, header = FALSE)

  # MGI Marker Accession ID | Marker Symbol | Status | Marker Name | cM Position | Chromosome | SWISS-PROT/TrEMBL Protein Accession IDs (space-delimited)
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

  # parse table from disk
  f = path_exists(f, NULL, try_compressed = TRUE)
  x = as_tibble(read_textfile_compressed(f, as_table = TRUE, datatable_skip = "GENE_RGD_ID"))
  # if(!file.exists(f)) {
  #   append_log(paste("file does not exist:", f), type = "error")
  # }
  # x = data.table::fread(f, data.table = FALSE, stringsAsFactors = FALSE, skip = "GENE_RGD_ID")

  # check that all required columns are present + rename columns
  cols = list(
    rgd_id = c("GENE_RGD_ID"),
    symbol = c("SYMBOL"),
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



#' helper function to match a protein table to HGNC via gene symbols (and synonyms)
#'
#' @description typically called from `export_stats_genesummary`
#' @param proteins_long long-format data.frame with columns protein_id and symbol pairs (i.e. no semicolon-delimited values)
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
#' @param use_synonyms boolean; use synonyms? if FALSE, only uses hgnc_symbol
idmap_symbols = function(proteins_long, hgnc, use_synonyms = TRUE) {
  if(!is.data.frame(proteins_long) || nrow(proteins_long) == 0 || !"symbol" %in% colnames(proteins_long) ) {
    append_log("proteins_long parameter must be a data.frame with column 'symbol'", type = "error")
  }
  if(!is.data.frame(hgnc) || nrow(hgnc) == 0 || !all(c("hgnc_id", "hgnc_symbol", "type", "value") %in% colnames(hgnc)) ) {
    append_log("hgnc parameter must be a data.frame with columns 'hgnc_id', 'hgnc_symbol', 'type', 'value' as typically prepared using the hgnc_lookuptable() function", type = "error")
  }
  if(length(use_synonyms) != 1 || !all(use_synonyms %in% c(TRUE, FALSE))) {
    append_log("use_synonyms parameter must be TRUE or FALSE", type = "error")
  }

  # initialize result column, if missing
  if(!"hgnc_id" %in% colnames(proteins_long)) {
    proteins_long$hgnc_id = NA_character_
  }

  # the hgnc mapping table contains official gene symbols in a separate column (i.e. not available as a "type")
  hgnc_subset_symbol = hgnc %>% distinct(hgnc_symbol, .keep_all = TRUE) # subset for efficiency

  # rows that we need to match by official gene symbol
  rows = is.na(proteins_long$hgnc_id) &
    # double-check that we don't even attempt to match missing gene symbols
    !is.na(proteins_long$symbol) & nchar(proteins_long$symbol) >= 2
  proteins_long$hgnc_id[rows] = hgnc_subset_symbol$hgnc_id[match(toupper(proteins_long$symbol[rows]), toupper(hgnc_subset_symbol$hgnc_symbol))]

  # fallback matching by synonym, using data in the "type" column
  if(use_synonyms) {
    rows = is.na(proteins_long$hgnc_id) & # analogous to above, but here update the selection of rows
      !is.na(proteins_long$symbol) & nchar(proteins_long$symbol) >= 2
    if(any(rows)) {
      hgnc_subset_synonym = hgnc %>% filter(type == "synonym") # subset only synonyms
      if(nrow(hgnc_subset_synonym) > 0) {
        proteins_long$hgnc_id[rows] = hgnc_subset_synonym$hgnc_id[match(toupper(proteins_long$symbol[rows]), toupper(hgnc_subset_synonym$value))]
      }
    }
  }

  return(proteins_long)
}



#' helper function to match a protein table to HGNC via an intermediate mapping table (e.g. MGI/RGD)
#'
#' @description typically called from `export_stats_genesummary`
#' @param proteins_long long-format data.frame with columns uniprot_id and symbol pairs (i.e. no semicolon-delimited values)
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
#' @param idmap MGI or RGD lookup table from `mgi_lookuptable()` or `rgd_lookuptable()`
#' @param use_symbols boolean; use symbols to map from proteins_long to the idmap table? if FALSE, only matches by uniprot identifier
idmap_uniprotid = function(proteins_long, hgnc, idmap, use_symbols = TRUE) {
  if(!is.data.frame(proteins_long) || nrow(proteins_long) == 0 || !all(c("uniprot_id", "symbol") %in% colnames(proteins_long)) ) {
    append_log("proteins_long parameter must be a data.frame with columns 'uniprot_id' and 'symbol'", type = "error")
  }
  if(!is.data.frame(hgnc) || nrow(hgnc) == 0 || !all(c("hgnc_id", "hgnc_symbol", "type", "value") %in% colnames(hgnc)) ) {
    append_log("hgnc parameter must be a data.frame with columns 'hgnc_id', 'hgnc_symbol', 'type', 'value' as typically prepared using the hgnc_lookuptable() function", type = "error")
  }
  if(!is.data.frame(idmap) || nrow(idmap) == 0 || ncol(idmap) < 3 || !all(c("symbol", "uniprot_id") %in% colnames(idmap)[-1]) || !colnames(idmap)[1] %in% hgnc$type) {
    append_log("idmap parameter must be a data.frame with columns 'symbol', 'uniprot_id' and the first column should contain identifiers that are also in the hgnc table (e.g. first column name should be 'mgi_id' or 'rgd_id'). e.g. as obtained from the mgi_lookuptable() function", type = "error")
  }
  if(length(use_symbols) != 1 || !all(use_symbols %in% c(TRUE, FALSE))) {
    append_log("use_symbols parameter must be TRUE or FALSE", type = "error")
  }

  # initialize result column, if missing
  if(!"hgnc_id" %in% colnames(proteins_long)) {
    proteins_long$hgnc_id = NA_character_
  }
  # overwrite xref column if exists
  proteins_long$xref_id = NA_character_

  # XREF ID @ first column of idmap
  XREF_ID_COLNAME = colnames(idmap)[1]
  idmap$XREF_ID_VALUES__TEMPCOL = unlist(idmap[,1], recursive = F, use.names = F) # ensure we're a vector and not a tibble of 1 column
  hgnc_subset_xref = hgnc %>% filter(type == XREF_ID_COLNAME)
  if(nrow(hgnc_subset_xref) == 0) {
    append_log(sprintf("assumed database type (from first idmap column): '%s'. However, values that match this crossreference database are not found as a 'type' in the provided hgnc table", colnames(idmap)[1]), type = "error")
  }


  ### first map from input table to idmap; prio matching by uniprot identifier
  rows = is.na(proteins_long$hgnc_id) & !is.na(proteins_long$uniprot_id)
  proteins_long$xref_id[rows] = idmap$XREF_ID_VALUES__TEMPCOL[match(proteins_long$uniprot_id[rows], idmap$uniprot_id)]

  # fallback to symbol matching (no `xref_id` yet)
  if(use_symbols) {
    rows = is.na(proteins_long$hgnc_id) & is.na(proteins_long$xref_id) & #
      # double-check that we don't even attempt to match missing gene symbols
      !is.na(proteins_long$symbol) & nchar(proteins_long$symbol) >= 2
    proteins_long$xref_id[rows] = idmap$XREF_ID_VALUES__TEMPCOL[match(toupper(proteins_long$symbol[rows]), toupper(idmap$symbol))]
  }


  ### map to HGNC using crossreference identifiers
  rows = is.na(proteins_long$hgnc_id) & !is.na(proteins_long$xref_id) # can be all FALSE, no problem
  proteins_long$hgnc_id[rows] = hgnc_subset_xref$hgnc_id[match(proteins_long$xref_id[rows], hgnc_subset_xref$value)]

  idmap$XREF_ID_VALUES__TEMPCOL = NULL
  return(proteins_long)
}



#' add hgnc_symbol and entrez/ensembl IDs if present in hgnc table
#'
#' @param x a data.table with column hgnc_id
#' @param hgnc HGNC lookup table from `hgnc_lookuptable()`
hgnc_add_xrefs = function(x, hgnc) {
  # symbol
  i = match(x$hgnc_id, hgnc$hgnc_id)
  x$symbol = hgnc$hgnc_symbol[i]

  # remove optional return values if present
  x$ensembl_id = NULL
  x$entrez_id = NULL

  # subset HGNC lookup table by respective crossref IDs
  hgnc_entrez = hgnc %>% filter(type == "entrez_id")
  hgnc_ensembl = hgnc %>% filter(type == "ensembl_id")

  if(nrow(hgnc_ensembl) > 0) {
    i = match(x$hgnc_id, hgnc_ensembl$hgnc_id)
    x$ensembl_id = hgnc_ensembl$value[i]
  }
  if(nrow(hgnc_entrez) > 0) {
    i = match(x$hgnc_id, hgnc_entrez$hgnc_id)
    x$entrez_id = hgnc_entrez$value[i]
  }

  return(x)
}
