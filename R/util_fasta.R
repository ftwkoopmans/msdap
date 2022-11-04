
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
#' @param x array of fasta headers
#' @param fasta_id_type fasta type. unused argument atm
fasta_header_to_gene = function(x, fasta_id_type = "uniprot") {
  res = sub(".* GN=([^ ]+).*", "\\1", x)
  res[res == x] = ""
  toupper(res)
}



#' Read non-redundant headers from all fasta files
#
#' @param fasta_files array of fill paths to fasta files
#' @param fasta_id_type fasta type, typically this is 'uniprot'
fasta_parse_file = function(fasta_files, fasta_id_type = "uniprot") {
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

  # distinct() at the end ensures we use data from the first fasta file for each unique accession, in case it's found in multiple
  tibble(
    idlong = fasta_header_to_id(fasta),
    gene = fasta_header_to_gene(fasta, fasta_id_type),
    header = fasta
  ) %>%
    mutate(idshort = fasta_id_short(idlong, fasta_id_type)) %>%
    distinct(idshort, .keep_all = T)
}
