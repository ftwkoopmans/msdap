
#' extract the protein identifier from a fasta header
#'
#' the first set of non-whitespace characters on the row is assumed to be the protein ID.
#' the leading '>' is not included (if present).
#'
#' this regex should be robust for all sorts of input, including those not following official standards.
#' should be ">proteinid whatever". But if additional whitespace is provided at the start we are robust to this.
#'
#' e.g. these should all yield 'pid' as the protein ID;
#' fasta_header_to_id(c(" > pid description", "> pid description", " pid description", "pid description"))
#'
#' @param x array of fasta headers
fasta_header_to_id = function(x) {
  sub("^[> ]*(\\S+).*", "\\1", x)
}



#' extract the "short ID" from a fasta header.
#'
#' for uniprot, a "short ID" is the accession between the pipe symbols.
#' for all other types, remove `<dot><record's version number>`.
#'
#' @param id array of fasta identifiers. eg; results from fasta_header_to_id()
#' @param fasta_id_type fasta type. if uniprot, the 'short ID' is the accession between the pipe symbols
fasta_id_short = function(id, fasta_id_type = "uniprot") {
  # for uniprot, a 'short ID' is the accession between the pipe symbols
  if(fasta_id_type == "uniprot") {
    return(sub(".*\\|(.*)\\|.*", "\\1", id))
  }
  # for all other types, remove <dot><record's version number> (only supports refseq as far as I know, but doesn't matter. should work for ensembl too because idshort==idlong )
  return(sub("\\.\\d*$", "", id))
}



#' extract gene symbol from fasta header string
#'
#' GN tag contains the gene symbol.
#' example; >sp|Q9Y2S6|TMA7_HUMAN Translation machinery-associated protein 7 OS=Homo sapiens (Human) OX=9606 GN=TMA7 PE=1 SV=1
#'
#' uniprot headers may contain "GN=-" or no GN tag at all !
#' example; >sp|Q6ZSR9|YJ005_HUMAN Uncharacterized protein FLJ45252 OS=Homo sapiens (Human) OX=9606 GN=- PE=2 SV=2
#' example; >tr|A2ALT2|A2ALT2_MOUSE Isoform of Q03288, Nonagouti (Fragment) OS=Mus musculus OX=10090 GN=a PE=4 SV=1
#' example; >sp|P15252|REF_HEVBR Rubber elongation factor protein OS=Hevea brasiliensis PE=1 SV=2
#'
#' @param x array of fasta headers
#' @param fasta_id_type fasta type. unused argument atm
fasta_header_to_gene = function(x, fasta_id_type = "uniprot") {
  # if multiple GN=XXX tags are present (should never occur for uniprot), the last element is taken
  res = sub(".* GN=([^ ]{2,}).*", "\\1", x) # require at least 2 characters
  res[res == x] = ""
  return(res)
}



#' Read non-redundant headers from all fasta files
#'
#' @description note that this function replaces semicolons in fasta headers with colons (since we use semicolons as delimiters upstream)
#' @param fasta_files array of fill paths to fasta files
#' @param fasta_id_type fasta type, typically this is 'uniprot'
#' @param uppercase_symbols convert all gene symbols to upper case? default: TRUE
fasta_parse_file = function(fasta_files, fasta_id_type = "uniprot", uppercase_symbols = TRUE) {
  # first, validate input files
  ff = NULL
  for (f in fasta_files) {
    f = path_exists(f, NULL, try_compressed = TRUE) # throws error if not found
    ff = c(ff, f)
  }

  # if we didn't abort with errors, we can start reading the validated file paths
  fasta = NULL
  for (f in ff) {
    l = read_textfile_compressed(f, as_table = F)
    char1 = substr(l, 1, 1)
    fasta = c(fasta, l[char1 == ">"]) # only fasta header lines
  }

  ### uniprot now includes semicolons in protein descriptions
  # importantly, this conflicts with using semicolons as a separation character
  # solution: replace all semicolons in fasta headers with colon symbols
  # example: >tr|B1AR09|B1AR09_MOUSE Isoform of B1AR10, Myeloid/lymphoid or mixed-lineage leukemia; translocated to, 6 OS=Mus musculus OX=10090 GN=Mllt6 PE=1 SV=1
  fasta = gsub(";", ":", fasta, fixed = TRUE)

  result = tibble(
    idlong = fasta_header_to_id(fasta),
    gene = fasta_header_to_gene(fasta, fasta_id_type),
    header = fasta
  ) %>%
    mutate(idshort = fasta_id_short(idlong, fasta_id_type)) %>%
    # distinct() at the end ensures we use data from the first fasta file for each unique accession, in case it's found in multiple
    # use the idlong because reverse/decoy tags might be prepended. e.g. "rev_tr|U3KQF4|U3KQF4_HUMAN"
    distinct(idlong, .keep_all = T) %>%
    mutate(isdecoy = fasta_identifier_isdecoy(idlong) | fasta_identifier_isdecoy(idshort))

  if(uppercase_symbols) {
    result$gene = toupper(result$gene)
  }
  return(result)
}



#' regex test if a protein identifier or fasta header is a decoy
#'
#' @examples \dontrun{
#' # decoy examples (should yield TRUE)
#' fasta_identifier_isdecoy(c(
#' ">rev_sp|A0A023PZB3|FMP49_YEAST Protein FMP49, mitochondrial ...",
#' "rev_sp|A0A023PYF4|YE145_YEAST",
#' "decoy_A0A023PZF2",
#' "P63044_decoy",
#' "P63044_rev",
#' "P63044_reverse",
#' ">sp|Q02956|KPCZ_MOUSE_decoy",
#' ">decoy_sp|Q02956|KPCZ_MOUSE"
#' ))
#' # non-decoy examples (should yield FALSE)
#' fasta_identifier_isdecoy(c(
#' ">sp|A0A023PZB3|FMP49_YEAST Protein FMP49, mitochondrial ...",
#' "sp|A0A023PYF4|YE145_YEAST",
#' "A0A023PZF2",
#' "P63044",
#' ">sp|Q02956|KPCZ_MOUSE"
#' ))
#' }
#' @param x strings to test
fasta_identifier_isdecoy = function(x) {
  !is.na(x) & is.character(x) & grepl("^>{0,1}REV_|_REV$|^>{0,1}REVERSE_|_REVERSE$|^>{0,1}DECOY_|_DECOY$", x, ignore.case = T)
}



#' lookup table for proteingroup identifiers (in long or short ID form) to idshort and isdecoy
#'
#' @examples \dontrun{
#' proteingroup_to_idshort_lookuptable(c(
#'   NA, "", # empty values are skipped
#'   "P36578", "P36578", "P36578", "P36578", # unique values are returned
#'   "P55011-3;P55011",
#'   "Q5TF21;rev_E9PJP2", "decoy_Q5TF21;E9PJP2",
#'   "sp|P63027|VAMP2_HUMAN;tr|K7ENK9|K7ENK9_HUMAN",
#'   "sp|P63027|VAMP2_HUMAN;rev_tr|K7ENK9|K7ENK9_HUMAN",
#'   "rev_sp|P63027|VAMP2_HUMAN;tr|K7ENK9|K7ENK9_HUMAN"
#' ))
#' }
#' @param x array of proteingroup identifiers
proteingroup_to_idshort_lookuptable = function(x) {
  x = unique(x)
  x = x[!is.na(x) & is.character(x)]
  if(length(x) == 0) return(NULL)

  tib_uprot = tibble(protein_id = x) %>%
    mutate(idlong = strsplit(protein_id, "[ ,;]+")) %>%
    # only empty strings are removed, which we don't have; tibble(id = letters[1:3], test = c("", "x;y", "z")) %>% mutate(test_list = strsplit(test, ";")) %>% unnest(cols = test_list, keep_empty = F)
    unchop(cols = idlong, keep_empty = FALSE)  %>%
    mutate(
      idshort = fasta_id_short(idlong),
      isdecoy = fasta_identifier_isdecoy(idlong) | fasta_identifier_isdecoy(idshort)
    )

  # rename decoy identifiers; strip current decoy flag (if any within idshort) and replace with standardized
  tmp = tib_uprot$idshort[tib_uprot$isdecoy]
  if(length(tmp) > 0) {
    tmp = gsub("^>{0,1}REV_|_REV$|^>{0,1}REVERSE_|_REVERSE$|^>{0,1}DECOY_|_DECOY$", "", tmp, ignore.case = T)
    tib_uprot$idshort[tib_uprot$isdecoy] = paste0("decoy_", tmp)
  }

  tib_uprot %>%
    group_by(protein_id) %>%
    summarise(idshort = paste(idshort, collapse = ";"),
              isdecoy = any(isdecoy)) %>%
    ungroup()
}


