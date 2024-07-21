
#' Create gene-level summary tables for your DEA and differential detection results
#'
#' @description
#'
#' To expedite downstream analysis, this function maps the proteingroup results from your
#' differential expression analysis to Human gene identifiers and exports these as Excel tables that are ready for
#' use with GO tools. Optionally, one can also include results from differential detection analysis.
#' Multiple proteingroups that map to the same unique gene are merged and there are various options
#' for dealing with ambiguous proteingroups.
#'
#' The resulting output tables can be used as input with the GOAT online geneset analysis tool @ \url{https://ftwkoopmans.github.io/goat/}
#'
#' *combining differential expression and differential detection stats*
#'
#' If differential detection (DD) was performed on this dataset and the `diffdetect_zscore_threshold` parameter is not NA,
#' the DD z-scores (per contrast) are converted to p-values by `pnorm(abs(zscore), lower.tail = F)`,
#' assuming these are approximately normally distributed (!), and then multiple testing correction by FDR is applied.
#' However, this is just an approximation; the differential detection scores are a very simplistic approach and should
#' be considered inferior to DEA analyses. Double-check the zscore distributions and treat these results with care !
#' Plot your DD results with `plot_differential_detect(dataset, zscore_threshold = 6)`.
#' Stringent cutoffs (absolute diffdetect_zscore_threshold of at least 5~6, or higher) are encouraged when using differential detection results.
#'
#' If DD is enabled (`diffdetect_zscore_threshold` parameter is not NA), the effectsizes from DEA and DD are both standardized to make
#' their distributions somewhat comparable. Note that the downside is that the effectsizes returned by this function
#' no longer will be the exact effectsizes returned by the DEA model (those can be then found @ output column `effectsize_dea`).
#'
#' *dealing with overlapping gene symbols between proteingroups*
#'
#' Importantly, there may be multiple proteingroups (isoforms) in any proteomics dataset that match the same gene symbol.
#' This function will return the proteingroup/row with the lowest p-value for each unique gene (per contrast, per DEA algorithm).
#' Your selected option for `gene_ambiguity` will affect this;
#' if for example "leading_gene" is selected as an option then a "GRIA1;GRIA2" proteingroup with pvalue=0.01 will be
#' chosen as representative value for "GRIA1" over a "GRIA1" proteingroup with pvalue=0.5 (see further below parameter documentation)
#'
#'
#' @examples \dontrun{
#'
#' ## note that in this example, we have renamed the downloaded files
#' ## to include a timestamp; useful for tracking versions
#' # load the HGNC data table you previously downloaded
#' hgnc_table = hgnc_lookuptable("C:/DATA/hgnc_complete_set__2024-01-01.txt")
#'
#' # MGI Mouse database, only needed for Mouse datasets
#' #mgi_table = mgi_lookuptable("C:/DATA/MRK_SwissProt_TrEMBL__2024-01-01.rpt")
#' # RGD Rat database, only needed for Rat datasets
#' #rgd_table = rgd_lookuptable("C:/DATA/GENES_RAT__2024-01-01.rpt")
#'
#' export_stats_genesummary(
#'   dataset,
#'   # set NA to ignore differential detection (default), or define an absolute
#'   # zscore threshold (e.g. 6). ! when using diffdetect, first review the
#'   # zscore distributions using MS-DAP function: plot_differential_detect()
#'   diffdetect_zscore_threshold = NA,
#'   # options for dealing with proteingroups that map to multiple genes:
#'   # leading_gene, prio_specific, only_specific
#'   gene_ambiguity = "prio_specific", # prio_specific and only_specific are recommended
#'   # always need the HGNC data table
#'   hgnc = hgnc_table,
#'   # example: for Mouse datasets, provide the MGI database for more accurate ID mapping
#'   # xref = mgi_table,
#'   # example: for Rat datasets, provide the RGD database for more accurate ID mapping
#'   # xref = rgd_table,
#'   # write the files to the same directory as main MS-DAP results
#'   output_dir = dataset$output_dir
#' )
#' }
#' @param dataset to use this function, your dataset should meet these requirements:
#' \itemize{
#' \item dataset was searched against a uniprot FASTA in MaxQuant/Spectronaut/DIA-NN/etc. (i.e. MS-DAP can only extract gene symbol info from from uniprot.org FASTA databases)
#' \item `import_fasta()` has been applied to this dataset
#' \item differential expression analysis has been applied (e.g. typical `analysis_quickstart()` workflow, or custom code that applied `filter_dataset()` and `dea()`)
#' }
#' @param gene_ambiguity options: leading_gene, prio_specific, only_specific  (default: prio_specific).
#' ## example table from GOAT online website
#' this list describes various scenarios; gene symbol(s) that represent some proteingroup together with a description of presented information
#' \itemize{
#' \item GRIA1 = protein group maps to exactly 1 gene
#' \item GRIA2 = protein group maps to exactly 1 gene
#' \item GRIA1;GRIA2 = ambiguous, one might want to use only the first entry ('leading' gene)
#' \item GRIA3;GRIA4 = ambiguous, but this row contributes a new gene (GRIA3)
#' \item tr|A8K0K0|A8K0K0_HUMAN;GRIA2;GRIA3 = ambiguous, but the first entry has no gene symbol only an accession
#' }
#'
#' ## options and their result when applied to above examples
#' \itemize{
#' \item leading_gene = use only the leading/first gene per proteingroup. e.g. "GRIA1;GRIA2" is mapped to "GRIA1" and will overwrite a specific "GRIA1" proteingroup if (and only if) the former has a lower p-value.
#'     all rows in this example will be mapped to a gene ID (respectively, GRIA1, GRIA2, GRIA1, GRIA3, GRIA2).
#' \item prio_specific = discard ambiguous proteingroups if their leading/first gene overlaps with a specific proteingroup. e.g. "GRIA1;GRIA2" is mapped to "GRIA1" only if there is no specific/unambiguous "GRIA1" proteingroup.
#'     rows 3 and 5 in this example are discarded because there exist specific entries for GRIA1 and GRIA2 (respectively, GRIA1, GRIA2, -, GRIA3, -). This approach favors rows that are specific and supplements the results only with ambiguous rows that contribute new information (leading gene is not in the results yet).
#' \item only_specific = all ambiguous proteingroups are discarded (included in output table, but human gene ID is set to NA).
#'     only the first 2 rows in this example are mapped (respectively, GRIA1, GRIA2, -, -, -).
#' }
#' @param diffdetect_zscore_threshold Optionally, a differential detection absolute z-score cutoff. Set to NA to disable diffdetect (default) to return only DEA results. When using this option, we strongly recommend to first review the zscore distributions using MS-DAP function: `plot_differential_detect()` -> double-check that your selected absolute zscore threshold cuts at intended positions in the distribution !
#' @param diffdetect_type type of differential detect scores. options:
#' 'auto' = set to 'detect' if this score is available, 'quant' otherwise
#' 'detect' = differential detection z-scores computed from only "detected" peptides (no MBR)
#' 'quant' = differential detection z-scores computed from all quantified peptides (uses MBR)
#' @param dea_logfc_instead_of_effectsize boolean parameter: should log2fc values be returned in the effectsize column? In some edge cases this might be helpful for downstream analysis but should generally be left at default value (FALSE)
#' @param output_dir where to write the Excel output tables (e.g. "C:/temp"). Set to NA to disable writing tables
#' @param hgnc HGNC lookup table, the output from function `hgnc_lookuptable()`
#' @param xref For Mouse and Rat datasets we recommend to provide MGI/RGD database files to increase ID mapping accuracy. For Mouse species, this should be the MGI lookup table from function `mgi_lookuptable()`. For Rat, the RGD lookup table from function `rgd_lookuptable()`
#' @param remove_nohgnc set to TRUE to remove all rows that could not be matched to a HGNC ID. default; FALSE
#' @returns the table(s) produced by this function are gene-level summaries of your statistical analyses.
#' The R data.table returned by this function is in long-format (i.e. results from multiple DEA algorithms and contrast are appended).
#' The Excel output files are split by contrast and DEA algorithm for convenience: these can be directly used as input for gene set analysis tools.
#' Proteingroups that failed to map to human genes are included in the output table but will have NA values in the hgnc_id column.
#'
#' A description of columns:
#' \itemize{
#' \item 'protein_id' contains the proteingroup identifier (concatenation of protein identifiers)
#' \item 'contrast' describes the statistical contrast (e.g. WT vs KO)
#' \item 'dea_algorithm' describes the DEA algorithm that configured. Note that if you select multiple algorithms while running `analysis_quickstart()`, there will be multiple rows in this table for each gene*contrast !
#' \item 'score_type' is either 'dea', 'dd' or 'dea,dd' indicating whether the effectsize/pvalue/pvalue_adjust shown for this row originates from differential expression analysis (DEA) or differential detection (DD), or both (i.e. statistics were available for both, in which case the metric with lowest pvalue was used to populate the effectsize/pvalue/pvalue_adjust columns)
#' \item 'peptide_count' describes the number of peptides that was used for DEA. So this is NOT the total number of detected peptides, rather this is the subset of peptides that passed your specified filtering rules (and were thus used for DEA)
#' \item 'log2fc' log2 foldchanges; for DEA results this value is taken straight from the DEA result table. For DD results this is the 'diffdetect zscore'
#' \item 'pvalue' analogous to `log2fc`
#' \item 'pvalue_adjust' analogous to `log2fc`
#' \item 'signif' boolean flag indicating significance. For DEA, the criteria for significance that were previously configured when performing DEA using `dea()` or
#' `analysis_quickstart()` functions are reused here (i.e. cutoffs for adjusted p-value and log2 foldchange). For DD results, any included protein (that passes the `diffdetect_zscore_threshold` cutoff) is set to TRUE.
#' \item The 'effectsize' column data depends on the user parameters/configuration;
#' When only differential detect is used, this contains the z-scores. When only DEA is used, this contains effectsizes from DEA as-is.
#' However, when DEA and differential detect are combined, the effectsizes from DEA are standardized (divided by std) such that effectsizes from
#' both statistics can be integrated into 1 distribution (and thus these can be compared/ranked in downstream analyses).
#' For rows that have both DD and DEA results, see item `score_type`
#' Not a perfect solution, but in combination with stringent `diffdetect_zscore_threshold` parameter this seems to work ok in our hands.
#' \item 'effectsize_dea' is only included when DD is enabled; it shows the DEA effectsizes as-is without any rescaling (and rescaled values merged with rescaled DD values are shown in 'effectsize')
#' \item 'gene_symbols_or_id' contains the uniprot gene symbols for respective accessions in the proteingroup (protein_id column) as per usual in MS-DAP.
#' \item 'hgnc_id' the mapped Human gene ID (can be NA in case of failed mappings)
#' \item 'symbol' the gene symbol of the respective 'hgnc_id' (can be NA in case of failed mappings)
#' \item 'ensembl_id' the Ensembl gene ID of the respective 'hgnc_id' (can be NA in case of failed mappings)
#' \item 'entrez_id' the NCBI Entrez gene ID of the respective 'hgnc_id' (can be NA in case of failed mappings)
#' \item 'gene' same as the 'entrez_id' column. Provided for drag-and-drop compatibility with the GOAT online tool
#' }
#' @export
export_stats_genesummary = function(dataset, gene_ambiguity = "prio_specific", diffdetect_zscore_threshold = NA, diffdetect_type = "auto", dea_logfc_instead_of_effectsize = FALSE, output_dir = NA, hgnc, xref = NULL, remove_nohgnc = FALSE) {
  result = NULL

  ### input validation
  if(!(is.list(dataset) && "de_proteins" %in% names(dataset))) {
    append_log("dataset parameter should be a dataset that contains DEA results (did you forget to setup contrasts and run the analysis_quickstart function, or dea() ?)", type = "error")
  }
  if(!(is.list(dataset) && "proteins" %in% names(dataset) && all(c("protein_id", "gene_symbols") %in% colnames(dataset$proteins)) && !all(dataset$proteins$protein_id == dataset$proteins$gene_symbols))) {
    append_log("dataset is missing essential information in the 'gene_symbols' column of the protein table. No FASTA has been imported for this dataset or the dataset was analyzed with an outdated MS-DAP version.\nTo use this function you won't have to re-run the entire analysis pipeline for this dataset, just run the import_fasta() function prior to applying export_stats_genesummary() to update the current dataset object (assuming this is a dataset searched against a uniprot FASTA)", type = "error")
  }
  if(!(length(gene_ambiguity) == 1 && all(gene_ambiguity %in% c("leading_gene", "prio_specific", "only_specific")))) {
    append_log("gene_ambiguity parameter should be 1 of the following options: leading_gene, prio_specific, only_specific", type = "error")
  }
  if(!(length(diffdetect_zscore_threshold) == 1 && (is.na(diffdetect_zscore_threshold) || (is.numeric(diffdetect_zscore_threshold) && is.finite(diffdetect_zscore_threshold) && diffdetect_zscore_threshold > 0)))) {
    append_log("diffdetect_zscore_threshold parameter should be either NA (to disable) or a positive number", type = "error")
  }
  if(!(length(dea_logfc_instead_of_effectsize) == 1 && dea_logfc_instead_of_effectsize %in% c(TRUE, FALSE))) {
    append_log("dea_logfc_instead_of_effectsize parameter should be TRUE or FALSE", type = "error")
  }
  if(!(length(output_dir) == 1 && (is.na(output_dir) || is.character(output_dir) && output_dir != ""))) {
    append_log("output_dir parameter should be either NA (to disable) or a (non-empty) character/string", type = "error")
  }
  if(!(is.data.frame(hgnc) && all(c("hgnc_id", "hgnc_symbol", "type", "value") %in% colnames(hgnc)))) {
    append_log("hgnc parameter should be a data.frame with columns hgnc_id, hgnc_symbol, type, value (output from hgnc_lookuptable() function)", type = "error")
  }
  if(!is.null(xref) && !(is.data.frame(xref) && nrow(xref) > 0 && ncol(xref) >= 3 && all(c("symbol", "uniprot_id") %in% colnames(xref)))) {
    append_log("xref parameter should be a data.frame with columns mgi_id/rgd_id, symbol, uniprot_id (output from mgi_lookuptable() or rgd_lookuptable() function)", type = "error")
  }
  if(!is.null(xref) && !colnames(xref)[1] %in% c("mgi_id", "rgd_id")) {
    append_log("xref parameter only works with MGI and RGD atm (output from mgi_lookuptable() or rgd_lookuptable() function): first column must be mgi_id or rgd_id", type = "error")
  }
  if(!(length(remove_nohgnc) == 1 && remove_nohgnc %in% c(TRUE, FALSE))) {
    append_log("remove_nohgnc parameter should be TRUE or FALSE", type = "error")
  }

  ### auto-detect species + issue user warning
  if("fasta_headers" %in% colnames(dataset$proteins)) {
    # example fasta header; ">tr|A0A024R4E5|A0A024R4E5_HUMAN Isoform of Q00341, High density lipoprotein binding protein (Vigilin)…"
    tmp = dataset$proteins %>% filter(protein_id %in% unique(dataset$de_proteins$protein_id)) %>% pull(fasta_headers)
    tmp = gsub(" .*", "", tmp) # remove description, full ID remains
    tmp = gsub(".*_", "", tmp) # remove everything except the species
    tmp = sort(table(tolower(tmp)), decreasing = TRUE) # counts for the unique species set
    tmp_subset = tmp[tmp > sum(tmp) * 0.6] # remove minor species

    # xref is provided: double-check that it matches species in the FASTA
    if(!is.null(xref)) {
      xref_is_mouse = colnames(xref)[1] == "mgi_id"
      xref_is_rat   = colnames(xref)[1] == "rgd_id"

      if(xref_is_mouse && ! "mouse" %in% names(tmp_subset)) {
        append_log(paste0("export_stats_genesummary: an MGI idmapping database is provided in parameter 'xref' while the vast majority of proteins in the dataset (auto-detected from uniprot headers in dataset$proteins) are not Mouse species.\nDetected majority species: ",
                          names(tmp_subset)[1], "\nOnly provide a crossreference database (parameter 'xref') when your dataset contains mostly proteins from that specific species!"), type = "error")
      }

      if(xref_is_rat && ! "rat" %in% names(tmp_subset)) {
        append_log(paste0("export_stats_genesummary: an RGD idmapping database is provided in parameter 'xref' while the vast majority of proteins in the dataset (auto-detected from uniprot headers in dataset$proteins) are not Rat species.\nDetected majority species: ",
                          names(tmp_subset)[1], "\nOnly provide a crossreference database (parameter 'xref') when your dataset contains mostly proteins from that specific species!"), type = "error")
      }
    } else {
      # no xref is provided, but maybe the user should !
      if(any(c("mouse", "rat") %in% names(tmp_subset))) {
        append_log(paste0("export_stats_genesummary: majority of protein species that was auto-detected from uniprot headers: ",
                          paste(unique(names(tmp_subset)), collapse = ","), "\nFor Mouse and Rat datasets we recommend to provide MGI/RGD database files to increase ID mapping accuracy (see further export_stats_genesummary() function documentation)"), type = "warning")
      }
    }
  }



  ### combine DEA and DD
  # Summarise the statistical data in your dataset, assuming DEA was applied.
  # note that if you analyzed multiple contrasts, optionally with multiple DEA algorithms, all results are appended in this table.
  allstats = summarise_stats(
    dataset,
    return_dea = TRUE,
    return_diffdetect = is.finite(diffdetect_zscore_threshold), # set to TRUE to integrate respective 'strong z-scores'
    dea_logfc_as_effectsize = dea_logfc_instead_of_effectsize,
    diffdetect_zscore_threshold = ifelse(is.finite(diffdetect_zscore_threshold), diffdetect_zscore_threshold, 6),
    diffdetect_type = diffdetect_type
  )



  ### create long-format table of uniprot_id*symbol for each unique proteingroup
  proteins = allstats %>%
    select(protein_id, gene_symbols) %>%
    distinct(protein_id, .keep_all = TRUE) %>%
    mutate(uniprot_id = strsplit(protein_id, "[ ,;]+"),
           symbol = strsplit(gene_symbols, "[ ,;]+"))

  if(!all(lengths(proteins$uniprot_id) == lengths(proteins$symbol))) {
    append_log("technical error: protein_id and gene_symbols columns contain different element counts. Rerun this dataset with the latest MS-DAP R function to amend, or create a GitHub ticket if this issue persists", type = "error")
  }

  proteins_long = proteins %>% select(-gene_symbols) %>% tidyr::unchop(cols = c(uniprot_id, symbol))



  ### map unique proteingroups to Human gene identifiers
  # first try to use the provided cross-reference database
  # (from uniprot to MGI, fallback from symbol to MGI, then MGI to HGNC)
  if(!is.null(xref)) {
    proteins_long = idmap_uniprotid(proteins_long, hgnc = hgnc, idmap = xref, use_symbols = TRUE)
  }
  # always map by gene symbol; will skip already-mapped (e.g. if mouse species & mapped via MGI)
  proteins_long = idmap_symbols(proteins_long, hgnc = hgnc, use_synonyms = TRUE)



  ### now join the unique proteingroup-to-HGNC mappings with table `allstats`
  # hgnc_list = a list column where each protein_id has an array of HGNC identifiers (length > 0), NA where missing/failed
  allstats = allstats %>%
    left_join(proteins_long %>% select(protein_id, hgnc_list = hgnc_id) %>% tidyr::chop(cols = "hgnc_list"), by = "protein_id")



  ### apply mapping rules
  ### importantly, sorting of `allstats` table prior to running below code MUST be ascending by pvalue

  allstats$hgnc_list_unique_notna = lapply(allstats$hgnc_list, function(g) { unique(na.omit(g)) }) # naive/slow, but prefer simple code
  # dev note: don't use `allstats %>% filter(lengths(hgnc_list) == 1)` because a proteingroup can contain multiple uniprot identifiers that all map to the same gene --> hgnc_list will not be length 1

  allstats$isambiguous = grepl(";", allstats$gene_symbols_or_id, fixed = TRUE) # infer ambiguous directly from uniprot gene symbols (i.e. failed mapping to HGNC have no effect)
  #  allstats$isambiguous = lengths(allstats$hgnc_list_unique_notna) > 1  # alternatively, determine "ambiguous" from successful HGNC mappings @ hgnc_list_unique_notna

  # take the first non-NA HGNC identifier (NA if none)
  allstats$hgnc_id = unlist(lapply(allstats$hgnc_list_unique_notna, function(g) {
    # always return exactly 1 value, thus downstream unlist() will yield an array that exactly matches table length
    if(length(g) == 0) return(NA_character_)
    return(g[1])
  }), recursive = FALSE, use.names = FALSE)

  if(gene_ambiguity == "only_specific") {
    # invalidate all entries that are ambiguous
    # importantly, do this before computing `hgnc_id_FORUNIQUE`
    allstats$hgnc_id[allstats$isambiguous] = NA_character_
  }

  if(gene_ambiguity == "prio_specific") {
    # sort table such that ambiguous hits are at the bottom; importantly, ties-of-ambiguity are still sorted as-is
    # taking first match per gene downstream = prio specific/unambiguous over ambiguous proteingroup that has a "better pvalue"
    allstats = allstats %>% arrange(isambiguous)
  }

  # ! importantly, add a unique ID to proteingroups not mapped to HGNC so we can use distinct() downstream without grouping NA values (failed mappings)
  allstats$hgnc_id_FORUNIQUE = ifelse(is.na(allstats$hgnc_id), allstats$gene_symbols_or_id, allstats$hgnc_id)
  # result table: retain unique rows per gene*contrast*dea_algorithm (assuming table has been sorted upstream)
  result = allstats %>% distinct(hgnc_id_FORUNIQUE, contrast, dea_algorithm, .keep_all = TRUE)



  ### finally, stable sort result table + add gene database crossreferences
  # sort table again ("prio_specific" option affects sorting, so here enforce stable ordering)
  result = result %>%
    arrange(desc(signif), pvalue) %>%
    select(-tidyselect::any_of(c("gene_symbols", "symbol", "hgnc_list", "hgnc_list_unique_notna", "hgnc_id_FORUNIQUE", "isambiguous"))) %>%
    hgnc_add_xrefs(hgnc = hgnc)
  # add named column specifically for GOAT
  result$gene = result$entrez_id



  tmp = allstats %>% distinct(protein_id, .keep_all = TRUE)
  log_pg_ambiguous = sum(tmp$isambiguous)

  if(gene_ambiguity == "only_specific") {
    # for "only_specific" option, don't count ambiguous genes @ success rate
    tmp = tmp %>% filter(isambiguous == FALSE)
    log_count_pg = nrow(tmp)
    log_count_pg_success = tmp %>% filter(!is.na(hgnc_id)) %>% nrow()

    log_string = sprintf(
      "%d / %d (%.1f%%) proteingroups were successfully mapped to %d unique human genes. %d additional ambiguous proteingroups were ignored (due to gene_ambiguity='only_specific')",
      log_count_pg_success, log_count_pg, log_count_pg_success / log_count_pg * 100, n_distinct(na.omit(tmp$hgnc_id)), log_pg_ambiguous
    )
  } else { # not 'only_specific'
    log_count_pg = nrow(tmp)
    log_count_pg_success = tmp %>% filter(!is.na(hgnc_id)) %>% nrow()
    log_string = sprintf(
      "%d / %d (%.1f%%) proteingroups were successfully mapped to %d unique human genes. ",
      log_count_pg_success, log_count_pg, log_count_pg_success / log_count_pg * 100, n_distinct(na.omit(tmp$hgnc_id))
    )

    tmp_g = tmp %>% filter(isambiguous == TRUE) %>% pull(gene_symbols_or_id)
    log_pg_ambiguous_overlap_unambiguous = sum(gsub(";.*", "", tmp_g) %in% tmp$gene_symbols_or_id)
    log_string = sprintf(
      "%s%d / %d ambiguous proteingroups had a leading gene symbol that overlapped with an unambiguous proteingroup",
      log_string, log_pg_ambiguous_overlap_unambiguous, log_pg_ambiguous
    )
  }
  append_log(log_string, type = "info")



  ### optionally, remove failed mappings
  if(remove_nohgnc) {
    result = result %>% filter(!is.na(hgnc_id))
  }



  ### print output tables
  if(is.na(output_dir)) {
    append_log("output_dir parameter set to NA, skip writing Excel output tables", type = "info")
  } else {
    for(contr in unique(result$contrast)) {
      for(algo_dea in unique(result$dea_algorithm)) {
        tmp = result %>% filter(contrast == contr & dea_algorithm == algo_dea)
        f = sprintf("%s/%s__%s.xlsx", output_dir, filename_strip_illegal_characters(sub("contrast: *", "", contr), strict = TRUE), algo_dea)
        openxlsx::write.xlsx(tmp, f)
      }
    }
  }

  return(invisible(result))
}
