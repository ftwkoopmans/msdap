
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
                                                          "** Experimental feature **",
                                                          "To easily prioritize proteins-of-interest which do not have sufficient abundance values for statistical analysis, but that do have many more detects in one sample group than the other, we provide a simple score based on the count data from the preceding columns.",
                                                          "This is an intentionally simplified approach that only serves to rank proteins for further qualitative analysis.",
                                                          "First, adjust the peptide detects per sample by the the total number in the sample to account for systematic differences.",
                                                          "The zscore metric in diff_detect_zscore_contrast: <contrast>'; A='count_peptides_detected_within_group_A';  B='count_peptides_detected_within_group_A';  result = scale(log2(B+min(Bi)) - log2(A+min(Ai))).    (where scale = subtract the mean, then divide by standard deviation)",
                                                          "",
                                                          "For all differential detection counts, be careful of over-interpretation and keep differences in sample group size (#replicates) and the total amount of peptides identified in a sample in mind !",
                                                          "",
                                                          "",
                                                          "Please refer to the online documentation for further details.")), keepNA = FALSE, headerStyle = header_style)
    openxlsx::writeData(wb, "statistics", tib, keepNA = FALSE, headerStyle = header_style)
    openxlsx::saveWorkbook(wb, file = fname_stats)
    # openxlsx::write.xlsx(tib, sprintf("%s/differential_abundance_analysis.xlsx", output_dir))
  } else {
    append_log("there are no statistical results in this dataset, nothing to write to Excel report", type = "info")
  }
}




#' placeholder title
#' @param peptides todo
#' @param col todo
#' @param ref_metadata todo
#' @param ref_colnames todo
tibble_prop_to_wide_and_align = function(peptides, col, ref_metadata, ref_colnames) {
  x = peptides %>%
    dplyr::select(peptide_id, sample_id, val=!!col) %>%
    tidyr::pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = val)

  # take peptides and their protein metadata from the protein intensity table, then align
  x = left_join(ref_metadata, x, by = "peptide_id")

  # now force aligned order and subset columns in detect tibble if we must
  x = x[ , intersect(ref_colnames, colnames(x))]
  return(x)
}



#' Create an Excel document with all peptide abundances
#'
#' Multiple sheets indicate results from all filters applied. For large datasets this results in huge files, so only recommended for small datasets.
#'
#' @param dataset your dataset
#' @param filename full file path to output file
#' @param norm_algorithm normalization algorithm(s)
#' @param norm_modebetween_protein_eset normalise protein-level data matrix used for DEA by 'modebetween' to correct for 'imbalance' between conditions introduced in peptide-to-protein rollup. Note that some algo_de, such as msqrob or msempire, directly operate on the peptide-level data thus are unaffected by this setting. Only useful for statistical models that first roll-up to protein-level, such as eBayes. default:FALSE
#'
#' @importFrom openxlsx createWorkbook addWorksheet createStyle writeData conditionalFormatting int2col saveWorkbook
#' @export
write_peptide_abundance_matrix_to_file = function(dataset, filename, norm_algorithm = "", norm_modebetween_protein_eset = F) {
  if(!grepl("\\.xlsx$", filename)) {
    append_log(paste("filename must end with .xlsx  @ ", filename), type = "error")
  }
  remove_file_if_exists(filename)

  ############ create peptide*sample intensity matrix
  tibw_intensity = peptide_abundance_table(dataset$peptides, dataset$proteins, dataset$samples, intensity_column_name = "intensity", norm_algorithm = norm_algorithm)
  # tibw_intensity = tibw_intensity %>% slice(10:100) %>% select(-WT_1) # debug; ensure subsets are nicely aligned

  ref_colnames = colnames(tibw_intensity)
  ref_metadata = tibw_intensity %>% select(!!c("peptide_id", colnames(dataset$proteins)))


  ############ create peptide*sample detection matrix
  tibw_detect = dataset$peptides %>%
    dplyr::select(peptide_id, sample_id, detect) %>%
    mutate(detect = as.numeric(detect)) %>%
    tidyr::pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = detect)

  # take peptides and their protein metadata from the protein intensity table, then align the peptide detect table
  tibw_detect = left_join(tibw_intensity %>% select(!!c("peptide_id", colnames(dataset$proteins))), tibw_detect,  by = "peptide_id") # %>% replace(is.na(.), FALSE)

  # now force aligned order and subset columns in detect tibble if we must
  tibw_detect = tibw_detect[,colnames(tibw_intensity)]


  ############ intersect both to get a matrix with only intensities whereever there is a detect
  colnames_data = setdiff(colnames(tibw_intensity), c("peptide_id", colnames(dataset$proteins)))

  tibw_intensity_hasdetect = tibw_intensity
  tmp_mask = as.matrix(tibw_detect[,colnames_data])
  tmp_mask[is.na(tmp_mask) | tmp_mask == F] = NA
  tibw_intensity_hasdetect[,colnames_data] = tibw_intensity_hasdetect[,colnames_data] * tmp_mask


  ############ create stylish Excel sheet that shows peptide abundances AND highlights the MBR values
  wb = openxlsx::createWorkbook()

  header_style = openxlsx::createStyle(textDecoration = "Bold")
  highlight_style = openxlsx::createStyle(bgFill = '#e0e0eb')

  openxlsx::addWorksheet(wb, "legend")
  openxlsx::addWorksheet(wb, "intensity")
  openxlsx::addWorksheet(wb, "intensity_where_identified")
  openxlsx::addWorksheet(wb, "identification")

  openxlsx::writeData(wb, "legend", data.frame(Legend=c("These data tables contain all peptide abundance values in the dataset.",
                                                        "If a normalization algorithms was specified, it is applied to the entire data matrix on sheet 2.",
                                                        "",
                                                        "Importantly, because no filtering is applied to this abundance table (prior to normalization),",
                                                        "these values may not be the exact same as those used in the differential abundance analysis !",
                                                        "eg; comparing sample groups A and B in a t-test, one typically filters for peptides found in 3+ replicates in both groups and applies normalization to that subset of the data, prior to statistical testing.",
                                                        "",
                                                        "Sheet 2 contains all normalized peptide intensities, color-coded if there is no matching identification (so for DDA these abundances originate from match-between-runs as there was no MS/MS identification. For DIA, these have confidence score < x).",
                                                        "Sheet 3 contains the subset of Sheet 2 for only those peptides that were confidently identified.",
                                                        "Sheet 4 shows a true/false flag indicating whether a peptide was confidently identified.",
                                                        "Following sheets contain any contrast-specific filtering+normalization, if user set parameter 'filter_within_contrast=TRUE'.",
                                                        "Finally, to make estimated protein abundance levels accessible we also include the data from sheets 2 and 3 rolled up to protein-level (summon all intensities per protein*sample).",
                                                        "",
                                                        "Please refer to the online documentation for further details.")), keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "intensity", tibw_intensity, keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "intensity_where_identified", tibw_intensity_hasdetect, keepNA = FALSE, headerStyle = header_style)
  openxlsx::writeData(wb, "identification", tibw_detect, keepNA = FALSE, headerStyle = header_style)

  ## conditional formatting approach
  # so it seems like openxlsx has issues creating across-sheet styles with conditional formatting
  # here we create a style, then manually update to inject the correct Excel compatible formula
  openxlsx::conditionalFormatting(wb, sheet = "intensity",
                                  cols = which(colnames(tibw_intensity) %in% colnames_data),
                                  rows = 2:(nrow(tibw_intensity)+1),
                                  # rule = 'B2="FALSE"',
                                  rule = 'AND(NOT(ISBLANK(###)),###=0)',
                                  style = highlight_style)

  # my workaround for across-sheet conditional formatting
  x = wb$worksheets[[2]]$conditionalFormatting
  # y = sub("B2=", paste0("identification!", openxlsx::int2col(which(colnames(tibw_intensity) %in% colnames_data)[1]), "2="), as.character(x), fixed = T)
  y = gsub("###", paste0("identification!", openxlsx::int2col(which(colnames(tibw_intensity) %in% colnames_data)[1]), "2"), as.character(x), fixed = T)
  names(y) = names(x)
  wb$worksheets[[2]]$conditionalFormatting <- y



  ############ intensities from user's filtering, as-is
  intensity_columns = grep("^intensity_", colnames(dataset$peptides), value=T)
  intensity_columns_shortname = sub("^intensity_contrast:\\s+", "", intensity_columns)
  is_toolong = nchar(intensity_columns_shortname) > 31
  intensity_columns_shortname[is_toolong] = paste0(substr(intensity_columns_shortname[is_toolong], 1, 26), " #", which(is_toolong))
  for(i in seq_along(intensity_columns)) {
    x = tibble_prop_to_wide_and_align(dataset$peptides, col=intensity_columns[i], ref_metadata = ref_metadata, ref_colnames=ref_colnames)
    openxlsx::addWorksheet(wb, intensity_columns_shortname[i])
    openxlsx::writeData(wb, intensity_columns_shortname[i], x, keepNA = FALSE, headerStyle = header_style)

    # tibw_int = peptides %>%
    #   dplyr::select(peptide_id, sample_id, intensity=!!(intensity_columns[i])) %>%
    #   tidyr::pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity)
    # # take peptides and their protein metadata from the protein intensity table, then align
    # tibw_int = left_join(tibw_intensity %>% select(!!c("peptide_id", colnames(proteins))), tibw_int,  by = "peptide_id")
    # # now force aligned order and subset columns in detect tibble if we must
    # tibw_int = tibw_int[,intersect(colnames(tibw_intensity), colnames(tibw_int))]
  }


  as_protein_eset = function(dataset, col_intensity, norm_modebetween_protein_eset) {
    eset_peptides = tibble_as_eset(dataset$peptides %>%
                                     select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, intensity=!!as.character(col_intensity)) %>%
                                     filter(is.finite(intensity)),
                                   dataset$proteins,
                                   dataset$samples)
    eset_proteins = eset_from_peptides_to_proteins(eset_peptides)

    if(norm_modebetween_protein_eset) {
      Biobase::exprs(eset_proteins) = normalize_matrix(Biobase::exprs(eset_proteins), algorithm = "modebetween", mask_sample_groups = Biobase::pData(eset_proteins)$group)
    }

    return(eset_proteins)
  }


  ############
  eset_peptides = tibble_as_eset(dataset$peptides %>% filter(is.finite(intensity) & detect == TRUE),
                                 dataset$proteins,
                                 dataset$samples)
  Biobase::exprs(eset_peptides) = normalize_matrix(Biobase::exprs(eset_peptides), algorithm = norm_algorithm, mask_sample_groups = MSnbase::pData(eset_peptides)$group)
  eset_proteins = eset_from_peptides_to_proteins(eset_peptides)
  if(norm_modebetween_protein_eset) {
    Biobase::exprs(eset_proteins) = normalize_matrix(Biobase::exprs(eset_proteins), algorithm = "modebetween", mask_sample_groups = Biobase::pData(eset_proteins)$group)
  }

  # take peptides and their protein metadata from the protein intensity table, then align
  tib = left_join(as_tibble(Biobase::exprs(eset_proteins)) %>% add_column(protein_id=rownames(Biobase::exprs(eset_proteins))), ref_metadata %>% select(-peptide_id) %>% distinct(protein_id, .keep_all = T), by = "protein_id") %>% arrange(gene_symbols_or_id)
  # now force aligned order and subset columns in detect tibble if we must
  tib = tib[ , intersect(ref_colnames, colnames(tib))]

  openxlsx::addWorksheet(wb, "protein_intensity_only_identify")
  openxlsx::writeData(wb, "protein_intensity_only_identify", tib, keepNA = FALSE, headerStyle = header_style)



  ############ protein intensities from 'all group filter'  and  'all intensities where detect'
  if("intensity_by_group" %in% colnames(dataset$peptides)) {
    eset_proteins = as_protein_eset(dataset, "intensity_by_group", norm_modebetween_protein_eset)
    # take peptides and their protein metadata from the protein intensity table, then align
    tib = left_join(as_tibble(Biobase::exprs(eset_proteins)) %>% add_column(protein_id=rownames(Biobase::exprs(eset_proteins))), ref_metadata %>% select(-peptide_id) %>% distinct(protein_id, .keep_all = T), by = "protein_id") %>% arrange(gene_symbols_or_id)
    # now force aligned order and subset columns in detect tibble if we must
    tib = tib[ , intersect(ref_colnames, colnames(tib))]

    ### v1
    # x = dataset$peptides %>%
    #   dplyr::select(peptide_id, protein_id, sample_id, val=intensity_by_group) %>%
    #   filter(!is.na(val)) %>%
    #   tidyr::pivot_wider(id_cols = c(peptide_id, protein_id), names_from = sample_id, values_from = val)
    # z = x %>% select(-peptide_id) %>% replace(is.na(.), 0) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA)
    # # y = as.data.frame(x %>% select(-peptide_id, -protein_id))
    # # Z = aggregate(y, by = list(protein_id=x$protein_id), FUN=sum, na.rm=T)
    # # all(z == Z, na.rm = T)
    #
    # # take peptides and their protein metadata from the protein intensity table, then align
    # tib = left_join(as_tibble(z) %>% replace(. == 0, NA), ref_metadata %>% select(-peptide_id) %>% distinct(protein_id, .keep_all = T), by = "protein_id") %>% arrange(gene_symbols_or_id)
    # # now force aligned order and subset columns in detect tibble if we must
    # tib = tib[ , intersect(ref_colnames, colnames(tib))]

    openxlsx::addWorksheet(wb, "protein_intensity_by_group")
    openxlsx::writeData(wb, "protein_intensity_by_group", tib, keepNA = FALSE, headerStyle = header_style)
  }

  if("intensity_all_group" %in% colnames(dataset$peptides)) {
    eset_proteins = as_protein_eset(dataset, "intensity_all_group", norm_modebetween_protein_eset)
    # take peptides and their protein metadata from the protein intensity table, then align
    tib = left_join(as_tibble(Biobase::exprs(eset_proteins)) %>% add_column(protein_id=rownames(Biobase::exprs(eset_proteins))), ref_metadata %>% select(-peptide_id) %>% distinct(protein_id, .keep_all = T), by = "protein_id") %>% arrange(gene_symbols_or_id)
    # now force aligned order and subset columns in detect tibble if we must
    tib = tib[ , intersect(ref_colnames, colnames(tib))]

    ### v1
    # x = dataset$peptides %>%
    #   dplyr::select(peptide_id, protein_id, sample_id, val=intensity_all_group) %>%
    #   filter(!is.na(val)) %>%
    #   tidyr::pivot_wider(id_cols = c(peptide_id, protein_id), names_from = sample_id, values_from = val)
    # z = x %>% select(-peptide_id) %>% replace(is.na(.), 0) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA)
    # # y = as.data.frame(x %>% select(-peptide_id, -protein_id))
    # # Z = aggregate(y, by = list(protein_id=x$protein_id), FUN=sum, na.rm=T)
    # # all(z == Z, na.rm = T)
    #
    # # take peptides and their protein metadata from the protein intensity table, then align
    # tib = left_join(as_tibble(z) %>% replace(. == 0, NA), ref_metadata %>% select(-peptide_id) %>% distinct(protein_id, .keep_all = T), by = "protein_id") %>% arrange(gene_symbols_or_id)
    # # now force aligned order and subset columns in detect tibble if we must
    # tib = tib[ , intersect(ref_colnames, colnames(tib))]

    openxlsx::addWorksheet(wb, "protein_intensity_all_group")
    openxlsx::writeData(wb, "protein_intensity_all_group", tib, keepNA = FALSE, headerStyle = header_style)
  }



  ## finally, write to file
  openxlsx::saveWorkbook(wb, file = filename, overwrite = TRUE)


  ## simplest implementation: write to Excel without any formatting
  # require(openxlsx)
  # list_of_datasets = list("peptide intensities" = tibw_intensity,
  #                         "peptide intensities where identified" = tibw_intensity_hasdetect,
  #                         "peptide identification" = tibw_detect)
  # write.xlsx(list_of_datasets, file = "C:/temp/test.xlsx")


  ## some reference code from earlier implementation:
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



#' placeholder title
#' for reference, standard pca (does not cope with missing values); summary(stats::prcomp(t(matrix_sample_intensities), center = T, scale. = F))
#' @param peptides todo
#' @param proteins todo
#' @param samples todo
#' @param intensity_column_name todo
#' @param norm_algorithm todo
peptide_abundance_table = function(peptides, proteins, samples, intensity_column_name = "intensity", norm_algorithm = "") {
  ## adopted from msdap::filter_dataset_by_group()
  mat = as_matrix_except_first_column(peptides %>%
                                        dplyr::select(peptide_id, sample_id, intensity = !!intensity_column_name) %>%
                                        dplyr::filter(!is.na(intensity)) %>%
                                        tidyr::pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity))
  for(alg in norm_algorithm) {
    if(alg == "") break
    mat = normalize_matrix(x_as_log2 = mat, algorithm = alg, mask_sample_groups = samples$group[match(colnames(mat), samples$sample_id)])
  }
  tib_norm = matrix_to_long(mat, value_name = "intensity", column_name = "sample_id", row_name = "peptide_id")
  ## adopted from msdap::filter_dataset_by_group()

  x = tib_norm
  x$protein_id = peptides$protein_id[match(x$peptide_id, peptides$peptide_id)]

  x_wide = tidyr::pivot_wider(x, id_cols = c(peptide_id, protein_id), names_from = sample_id, values_from = intensity)
  x_wide = x_wide[ , order(match(colnames(x_wide), c("peptide_id", "protein_id", samples$sample_id)))]
  # boxplot(as.matrix(x_wide[,colnames(x_wide) %in% samples$sample_id]), outline=F, las=1, horizontal=T, cex.axis=.5)

  # add protein metadata and sort by {gene, peptide}
  right_join(proteins, x_wide, by = "protein_id") %>% arrange(gene_symbols_or_id, peptide_id)
}
