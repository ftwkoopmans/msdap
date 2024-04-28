
#' write statistical results to file
#'
#' combines both the DEA and "differential testing" results into output file differential_abundance_analysis.xlsx
#'
#' @param dataset your dataset
#' @param output_dir path to existing directory
#'
#' @importFrom openxlsx createWorkbook addWorksheet createStyle writeData saveWorkbook write.xlsx
#' @export
export_statistical_results = function(dataset, output_dir) {
  fname_stats = path_append_and_check(output_dir, "differential_abundance_analysis.xlsx")
  remove_file_if_exists(fname_stats)

  tib = tib_dea = dea_results_to_wide(dataset)
  tib_dt = diffdetect_results_to_wide(dataset)
  if(nrow(tib_dea) == 0) {
    tib = tib_dt
  }
  if(nrow(tib_dea) > 0 && nrow(tib_dt) > 0) {
    tib = full_join(tib_dea, tib_dt, by="protein_id")
  }

  if(nrow(tib) > 0) {
    # add # unique peptides per protein in entire dataset
    pepcount = dataset$peptides %>% select(peptide_id, protein_id) %>% distinct(peptide_id, .keep_all = T) %>% count(protein_id, name = "unique_peptides")
    tib = add_column(tib, unique_peptides = pepcount$unique_peptides[data.table::chmatch(tib$protein_id, pepcount$protein_id)], .before = 1)

    # if there is protein metadata, join and put those columns in front
    if(is_tibble(dataset$proteins) && nrow(dataset$proteins) > 0) {
      tib = full_join(tib, dataset$proteins, by="protein_id") %>% select(!!colnames(dataset$proteins), everything())
    }

    # write as Excel file
    wb = openxlsx::createWorkbook()
    header_style = openxlsx::createStyle(textDecoration = "Bold")
    openxlsx::addWorksheet(wb, "legend")
    openxlsx::addWorksheet(wb, "statistics")
    openxlsx::writeData(wb, "legend", data.frame(Legend=c("This data table contains the results from statistical analyses, each row describes one protein(group).",
                                                          "",
                                                          "The first few columns describe the protein metadata; the protein accessions together with the respective descriptions and gene symbol(s) as listed in the provided fasta file. If no gene symbol is provided in the fasta file for a protein, the 'gene_symbols_or_id' column simply contains the value from 'protein_id'.",
                                                          "",
                                                          "The next set of columns describe the results from each contrast (comparison of sample groups).",
                                                          "'peptides_used_for_dea' shows the number of peptides that pass the filter rules within some contrast (eg; comparing WT vs KO, how many peptides for protein P could be used for statistical analysis?).",
                                                          "",
                                                          "For each statistical model applied (DEqMS, MSqRob, etc.), the log2 foldchange, pvalue and qvalue (pvalue adjusted for multiple testing) are provided.",
                                                          'Note that a MS-DAP contrast for "A vs B" returns foldchanges for B/A. For example, for the contrast "control vs disease" a positive log2 foldchange implies protein abundances are higher in the "disease" sample group.',
                                                          "If multiple DEA models were applied, the column signif_count_<contrast> shows how many algorithms deemed each protein significant.",
                                                          "If you choose to post-hoc combine results from multiple DEA models, we recommend selecting proteins significant in at least 2 of these algorithms as a compromise between using multiple algorithms to maximize your search for proteins-of-interest and robustness against false positives (eg; proteins uniquely scoring well in one statistical model and not the others).",
                                                          "Note; if you apply multiple DEA algorithms but don't follow this selection rule (eg; protein found in N+ DEA algorithms), you should apply some form of multiple-testing-correction to compensate for running multiple hypotheses tests. Consult your local statistician.",
                                                          "",
                                                          "After the differential expression analysis results, qualitative data shows in how many replicates within each group peptides were identified (MS/MS for DDA, identification confidence < x for DIA).",
                                                          "'count_peptides_detected_within_group_<sample group>' shows how many time a peptide for the protein on this row was detected over all replicates within a group (for each sample, count unique detected peptides, then sum for all replicates within group). So if a protein has 3 peptides in total and there are 4 replicates in a group, this is at most 12.",
                                                          "This is pseudo-count data, intended to be used for qualitative analysis when peptide abundances for a protein are absent in one sample group (or quantified in too few replicates to use for differential expression analysis).",
                                                          "For instance, in a 'WT vs KO' comparison one should see no (or very few) such peptide detections in the KO for the target protein and many more in the WT sample group.",
                                                          "",
                                                          "Differential Detection:",
                                                          "Compute a z-score for each protein based on the total number of detected peptides per sample group.",
                                                          "To easily prioritize proteins-of-interest, some which may not have sufficient abundance values for differential expression analysis but that do have many more detects in one sample group than the other, we provide a simple score based on identification count data.",
                                                          "This is a simplified approach that is intended to rank proteins for further qualitative analysis, be careful of over-interpretation and keep differences in sample group size (#replicates) and the absolute amount of peptides identified in a sample in mind !",
                                                          "columns that start with 'diff_detect_zscore_contrast:' contain the z-scores for each contrast (only for proteins observed in at least N samples samples in either group, to exclude one-hit-wonders)",
                                                          "",
                                                          "Please refer to the online documentation for further details: https://github.com/ftwkoopmans/msdap")), keepNA = FALSE, headerStyle = header_style)

    openxlsx::writeData(wb, "statistics", tib, keepNA = FALSE, headerStyle = header_style)
    openxlsx::saveWorkbook(wb, file = fname_stats)
  } else {
    append_log("there are no statistical results in this dataset, nothing to write to Excel report", type = "info")
  }
}



#' Rollup and export to protein-level TSV table, applied to all peptide-level intensity data
#'
#' generated a protein-level abundance matrix for each column in dataset$peptides that starts with "intensity"
#'
#' @param dataset your dataset
#' @param rollup_algorithm strategy for combining peptides to proteins. For valid options, see rollup_pep2prot()
#' @param output_dir output directory where all output files are stored, must be an existing directory
#' @export
export_protein_abundance_matrix = function(dataset, rollup_algorithm, output_dir) {
  stopifnot(dir.exists(output_dir))
  # iterate all intensity data; input data as-is (intensity) and filtered+normalized data, if available, from all_group or by_contrast filters
  for(col_contr_intensity in sort(grep("^intensity", colnames(dataset$peptides), value = T))) {

    # protein-level data matrix
    m = rollup_pep2prot(dataset$peptides %>%
                          select(sample_id, protein_id, peptide_id, intensity=!!as.character(col_contr_intensity)) %>%
                          filter(is.finite(intensity)),
                        intensity_is_log2 = T, rollup_algorithm = rollup_algorithm, return_as_matrix = T)
    # sort columns by order in sample metadata table. can safely match because each sample_id in dataset$peptides must be present in dataset$samples (provided input dataset is 'valid')
    m = m[ , order(match(colnames(m), dataset$samples$sample_id)), drop=F]

    # add protein metadata
    tib = dataset$proteins %>% inner_join(as_tibble(m) %>% add_column(protein_id = rownames(m)), by="protein_id") %>% arrange(protein_id)
    if("accessions" %in% colnames(tib)) {
      tib = tib %>% select(-accessions) # not useful for user, redundant with protein_id column in virtually all datasets
    }
    if("gene_symbols_or_id" %in% colnames(tib)) {
      tib = tib %>% arrange(gene_symbols_or_id!=protein_id, gene_symbols_or_id) # proteins without gene symbol first, then sort by symbol
    }

    ## write to file
    # generate filename. if very long (eg; huge contrast name + long path in output_dir), try to shorting with md5 hash
    fname = path_append_and_check(output_dir, sprintf("protein_abundance__%s.tsv", column_intensity_to_label(col_contr_intensity)))
    # check if output file is writable
    remove_file_if_exists(fname)
    # write to TSV file. note that we enforce a dot as decimal character (so ignore user's locale)
    utils::write.table(tib, file = fname, quote = F, sep="\t", na="", row.names = F, col.names = T, dec = ".")
  }
}



#' Export peptide-level data to TSV files
#'
#' Generates separate files for input data as-is (i.e. all peptides, intensity values as provided in the imported dataset)
#' and all filtering&normalization applied to the dataset.
#'
#' i.e. 1 output file for every column in dataset$peptides that starts with "intensity"
#'
#' @param dataset your dataset
#' @param output_dir output directory where all output files are stored, must be an existing directory
#' @export
export_peptide_abundance_matrix = function(dataset, output_dir) {
  stopifnot(dir.exists(output_dir))
  # iterate all intensity data; input data as-is (intensity) and filtered+normalized data, if available, from all_group or by_contrast filters
  for(col_contr_intensity in sort(grep("^intensity", colnames(dataset$peptides), value = T))) {
    # peptide-level data matrix
    tib = dataset$peptides %>%
      select(sample_id, protein_id, peptide_id, intensity=!!as.character(col_contr_intensity)) %>%
      filter(is.finite(intensity))

    # samples with values in current intensity column, ordered like sample metadata table
    sid = intersect(dataset$samples$sample_id, unique(tib$sample_id))

    tibw = tib %>%
      pivot_wider(id_cols = c(protein_id, peptide_id), names_from = "sample_id", values_from = "intensity") %>%
      # add protein metadata
      left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by = "protein_id") %>%
      # sort columns by order in sample metadata table
      select(gene_symbols_or_id, protein_id, peptide_id, !!sid) %>%
      # sort output table
      arrange(gene_symbols_or_id!=protein_id, gene_symbols_or_id, protein_id, peptide_id) # proteins without gene symbol first, then sort by symbol

    ## write to file
    # generate filename. if very long (eg; huge contrast name + long path in output_dir), try to shorting with md5 hash
    fname = path_append_and_check(output_dir, sprintf("peptide_abundance__%s.tsv", column_intensity_to_label(col_contr_intensity)))
    # check if output file is writable
    remove_file_if_exists(fname)
    # write to TSV file. note that we enforce a dot as decimal character (so ignore user's locale)
    utils::write.table(tibw, file = fname, quote = F, sep="\t", na="", row.names = F, col.names = T, dec = ".")
  }

}
