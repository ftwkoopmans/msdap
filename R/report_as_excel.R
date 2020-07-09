
#' write_statistical_results_to_file
#'
#' @param dataset your dataset
#' @param output_dir path to existing directory
#'
#' @importFrom openxlsx createWorkbook addWorksheet createStyle writeData saveWorkbook write.xlsx
#' @export
write_statistical_results_to_file = function(dataset, output_dir) {
  fname_stats = sprintf("%s/differential_abundance_analysis.xlsx", output_dir)
  remove_file_if_exists(fname_stats)

  tib = tib_dea = dea_results_to_wide(dataset)
  tib_dt = diffdetect_results_to_wide(dataset)
  if(nrow(tib_dea) == 0) {
    tib = tib_dt
  }
  if(nrow(tib_dea) > 0 && nrow(tib_dt) > 0) {
    tib = full_join(dea_results_to_wide(dataset),
                    diffdetect_results_to_wide(dataset),
                    by="protein_id")
  }


  if(nrow(tib) > 0) {
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
                                                          "The first few columns describe the protein metadata; the protein accessions together with the respective descriptions and gene symbol(s) as listed in the provided fasta file. If no gene symbol is given in the fasta file for a protein, the 'gene_symbols_or_id' column simply contains the value from 'protein_id'.",
                                                          "",
                                                          "The next set of columns describe the results from each contrast (comparison of sample groups).",
                                                          "'peptides_used_for_dea' shows the number of peptides that pass the filter rules within some contrast (eg; comparing WT vs KO, how many peptides for protein P could be used for statistical analysis?).",
                                                          "",
                                                          "For each statistical model applied (eBayes, msqrob, msempire), the log2 foldchange, pvalue and qvalue (adjusted pvalue) are provided.",
                                                          "If multiple DEA models were applied, the column signif_count_<contrast> shows how many algorithms deemed each protein significant.",
                                                          "As shown in our data analyses, applying eBayes/msqrob/msempire and then selecting proteins significant in at least 2 of these algorithms is a compromise between using multiple algorithms to maximize your search for proteins-of-interest and robustness against false positives (eg; proteins uniquely scoring well in one statistical model and not the others).",
                                                          "Note; if you apply multiple DEA algorithms but don't follow this selection rule (eg; protein found in N+ DEA algorithms), you should apply some form of multiple-testing-correction to compensate for running multiple hypotheses tests. Consult your local statistician.",
                                                          "",
                                                          "After the differential expression analysis results, qualitative data shows in how many replicates within each group peptides were identified (MS/MS for DDA, identification confidence < x for DIA).",
                                                          "'count_peptides_detected_within_group_<sample group>' shows how many time a peptide for the protein on this row was detected over all replicates within a group (for each sample, count unique detected peptides, then sum for all replicates within group). So if a protein has 3 peptides in total and there are 4 replicates in a group, this is at most 12.",
                                                          "This is pseudo-count data, intended to be used for qualitative analysis when peptide abundances for this protein are absent in one sample group (or quantified in too few replicates to use for differential abundance analysis).",
                                                          "For instance, in a 'WT vs KO' comparison one should see no (or very few) such peptide detections in the KO for the target protein and many more in the WT sample group.",
                                                          "",
                                                          "Differential Detection:",
                                                          "Compute a z-score for each protein based on the total number of detected peptides per sample group.",
                                                          "To easily prioritize proteins-of-interest, some which may not have sufficient abundance values for differential expression analysis but that do have many more detects in one sample group than the other, we provide a simple score based on identification count data.",
                                                          "This is a simplified approach that is intended to rank proteins for further qualitative analysis, be careful of over-interpretation and keep differences in sample group size (#replicates) and the absolute amount of peptides identified in a sample in mind !",
                                                          "columns that start with 'diff_detect_zscore_contrast:' contain the z-scores for each contrast, proteins that meet pre-defined selection criteria are flagged in columns that start with 'diff_detect_zscore_candidate_contrast:'",
                                                          "",
                                                          "Please refer to the online documentation for further details: https://github.com/ftwkoopmans/msdap")), keepNA = FALSE, headerStyle = header_style)

    openxlsx::writeData(wb, "statistics", tib, keepNA = FALSE, headerStyle = header_style)
    openxlsx::saveWorkbook(wb, file = fname_stats)
  } else {
    append_log("there are no statistical results in this dataset, nothing to write to Excel report", type = "info")
  }
}



#' Create an Excel document with peptide and protein abundances
#'
#' Conditional formatting is used to highlight match-between-runs peptides (low confidence data points).
#'
#' Keep in mind that large datasets may result in (too) large Excel files, so in some cases you're better of pre-processing your data in R.
#'
#' included tables;
#' - peptide intensity (-filter +norm. conditional formatting)
#' - peptide identification (1/0 flags)
#' - peptide intensity_all_group (+filter +norm)
#' - protein roll-up from peptide intensity
#' - protein roll-up from peptide intensity_all_group
#'
#' @param dataset your dataset
#' @param filename full file path to output file
#' @param norm_algorithm normalization algorithm(s)
#'
#' @importFrom openxlsx createWorkbook addWorksheet createStyle writeData conditionalFormatting int2col saveWorkbook
#' @importFrom data.table chmatch
#' @export
write_peptide_abundance_matrix_to_file = function(dataset, filename, norm_algorithm = "") {
  if(!grepl("\\.xlsx$", filename)) {
    append_log(paste("filename must end with .xlsx  @ ", filename), type = "error")
  }
  remove_file_if_exists(filename)
  start_time <- Sys.time()


  ############ peptides

  # normalize the unfiltered peptide data
  if(all(norm_algorithm != "")) {
    dataset$peptides$intensity_norm = normalize_intensities(data.table::setDT(dataset$peptides %>% select(key_peptide_sample, key_peptide, key_protein, key_sample, key_group, intensity)), norm_algorithm = norm_algorithm)
  } else {
    dataset$peptides$intensity_norm = x$intensity
  }

  # from long-format peptide tibble to a sorted wide-format table that includes protein metadata (only protein_id and gene symbols, no fasta headers to save disk space at huge datasets)
  tbl_pep = dataset$peptides %>%
    dplyr::select(protein_id, peptide_id, sample_id, intensity_nofilter = intensity_norm, intensity_filter_allgroup = intensity_all_group, detected = detect) %>%
    # samples (columns) ordered like dataset$samples tibble
    dplyr::arrange(data.table::chmatch(sample_id, dataset$samples$sample_id)) %>%
    dplyr::mutate(detected = as.integer(detected)) %>%
    tidyr::pivot_wider(id_cols = c(protein_id, peptide_id), names_from = sample_id, values_from = c(intensity_nofilter, intensity_filter_allgroup, detected), names_sep = ": ") %>%
    dplyr::left_join(dataset$proteins %>% select(-accessions, -fasta_headers), by = "protein_id") %>%
    dplyr::select(protein_id, gene_symbols_or_id, tidyr::everything()) %>%
    # peptides (rows) ordered by gene symbol, then protein group, then peptide ID
    dplyr::arrange(gene_symbols_or_id, protein_id, peptide_id)


  ############ proteins

  # peptide to protein roll-up within dplyr
  tbl_prot = dataset$peptides %>%
    dplyr::select(sample_id, protein_id, intensity_norm, intensity_all_group) %>%
    dplyr::mutate(intensity_norm = 2^intensity_norm, intensity_all_group = 2^intensity_all_group) %>%
    dplyr::group_by(sample_id, protein_id) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate(intensity_norm = log2(intensity_norm), intensity_all_group = log2(intensity_all_group))

  # analogous to peptide table
  tbl_prot = tbl_prot %>%
    dplyr::select(sample_id, protein_id, intensity_nofilter = intensity_norm, intensity_filter_allgroup = intensity_all_group) %>%
    # samples (columns) ordered like dataset$samples tibble
    dplyr::arrange(data.table::chmatch(sample_id, dataset$samples$sample_id)) %>%
    tidyr::pivot_wider(id_cols = protein_id, names_from = sample_id, values_from = c(intensity_nofilter, intensity_filter_allgroup), names_sep = ": ") %>%
    dplyr::left_join(dataset$proteins, by = "protein_id") %>%
    dplyr::select(!!colnames(dataset$proteins), tidyr::everything()) %>%
    # proteins (rows) ordered by gene symbol, then protein group
    dplyr::arrange(gene_symbols_or_id, protein_id)


  ############ create stylish Excel sheet that shows peptide abundances AND highlights the MBR values

  wb = openxlsx::createWorkbook()

  header_style = openxlsx::createStyle(textDecoration = "Bold")
  highlight_style = openxlsx::createStyle(bgFill = '#e0e0eb')

  openxlsx::addWorksheet(wb, "legend")
  openxlsx::addWorksheet(wb, "peptides")
  openxlsx::addWorksheet(wb, "proteins")

  openxlsx::writeData(wb, "legend", data.frame(Legend=c("These data tables contain normalized peptide abundance values.",
                                                        "",
                                                        "Sheet 2 contains peptide-level data.",
                                                        "Columns that start with 'intensity_nofilter:' are log2 intensity values without any filtering (e.g. all data points available in the input dataset).",
                                                        "Columns that start with 'intensity_filter_allgroup:' are log2 intensity values for those peptides that pass the selection criteria in all sample groups. These procedures follow user settings, such as the 'min_detect' parameter that dictates in how many sample a peptide must be confidently identifier per sample group.",
                                                        "Columns that start with 'detected:' indicate peptides that were confidently detected with a 1 and 'low confidence' data points with a 0. So for DDA experiments such 'low confidence' abundance values originate from match-between-runs as there was no MS/MS identification. For DIA, these do not meet the confidence score threshold that was set while importing the dataset in MS-DAP.",
                                                        "Color-coding is used to highlight the 'low confidence' data points where detection/identification is lacking.",
                                                        "",
                                                        "Sheet 3 contains protein abundances values that were computed by roll-up from the data shown in the peptide sheet. Summation of all (non-log) peptide intensities per protein is used for peptide-to-protein roll-up.",
                                                        "",
                                                        "Note that if your DEA used the filtering and normalization 'by contrast' this data does not exactly match the input values for DEA (you can always find the exact data used for DEA in the long-format data table dataset$peptides).",
                                                        "",
                                                        "Please refer to the online documentation for further details: https://github.com/ftwkoopmans/msdap")), keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "peptides", tbl_pep, keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "proteins", tbl_prot, keepNA = FALSE, headerStyle = header_style)


  ## conditional formatting approach
  # so it seems like openxlsx has issues creating across-sheet styles with conditional formatting
  # here we create a style, then manually update to inject the correct Excel compatible formula
  openxlsx::conditionalFormatting(wb,
                                  sheet = "peptides",
                                  cols = grep("^intensity_nofilter", colnames(tbl_pep)),
                                  rows = 2:(nrow(tbl_pep)+1),
                                  rule = 'AND(NOT(ISBLANK(###)),###=0)',
                                  style = highlight_style)

  openxlsx::conditionalFormatting(wb,
                                  sheet = "peptides",
                                  cols = grep("^intensity_filter_allgroup", colnames(tbl_pep)),
                                  rows = 2:(nrow(tbl_pep)+1),
                                  rule = 'AND(NOT(ISBLANK(@@@)),***=0)',
                                  style = highlight_style)

  # my workaround for across-sheet conditional formatting
  x = wb$worksheets[[2]]$conditionalFormatting
  # y = gsub("###", paste0(openxlsx::int2col(grep("^detected", colnames(tbl_pep))[1]), "2"), as.character(x), fixed = T)
  y = gsub("###", paste0("peptides!", openxlsx::int2col(grep("^detected", colnames(tbl_pep))[1]), "2"), as.character(x), fixed = T)
  y = gsub("@@@", paste0("peptides!", openxlsx::int2col(grep("^intensity_filter_allgroup", colnames(tbl_pep))[1]), "2"), as.character(y), fixed = T)
  y = gsub("***", paste0("peptides!", openxlsx::int2col(grep("^detected", colnames(tbl_pep))[1]), "2"), as.character(y), fixed = T)
  names(y) = names(x)
  wb$worksheets[[2]]$conditionalFormatting <- y


  ############  finally, write to file
  openxlsx::saveWorkbook(wb, file = filename, overwrite = TRUE)

  append_log_timestamp("Excel document with peptide and protein abundances", start_time)


  ############ some reference code from earlier implementation:
  ## manually color-code each cell. advantage: full control within R. disadvantage: slow to write workbook to file
  # mbrstyle = openxlsx::createStyle(fgFill = "#e0e0eb") # https://github.com/awalker89/openxlsx/issues/197
  # df_int = as.matrix(tibw_intensity[,colnames_data])
  # df_det = as.matrix(tibw_detect[,colnames_data])
  # index = which(!is.na(df_int) & !df_det, arr.ind = T)
  # column_offset_metadata = which(colnames(tibw_intensity) %in% colnames_data)[1] -1
  # for(row in 1:nrow(index)) {
  #   i = index[row,1] + 1
  #   j = index[row,2] + column_offset_metadata
  #   openxlsx::addStyle(wb, sheet = 1, style = mbrstyle, rows = i, cols = j, gridExpand = FALSE)
  # }
}
