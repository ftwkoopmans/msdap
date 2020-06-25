
prepare_test_datasets = function(subset_n_rows = NA) {
  ### example dataset and application of MS-EmpiRe following the documentation at: https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R
  f <- system.file("extdata", "c1_c3.data", package = "msEmpiRe")
  p <- system.file("extdata", "c1_c3.pdata", package = "msEmpiRe")
  # loading data from installed data sets. The dataset is the yeast spike-in benchmarking data of O'Connell et al., as discussed in the paper
  suppressMessages(data <- msEmpiRe::read.standard(f, p,
                                                   prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                                   signal_pattern="c.*rep.*"))
  # extract the first two conditions
  conditions <- msEmpiRe::extract_conditions(data)
  conditions <- conditions[, c(1,2)]
  # removing peptides that are detected in less than 2 samples per condition
  tmp = capture.output(msempire_eset_filtered <- msEmpiRe::filter_detection_rate(data, condition=conditions))
  # optionally, select the first N peptides (that passed basic filter criteria)
  if(!is.na(subset_n_rows)) {
    msempire_eset_filtered = msempire_eset_filtered[1:subset_n_rows, ]
  }

  tmp = log2(Biobase::exprs(msempire_eset_filtered))
  tmp[!is.finite(tmp)] = NA
  Biobase::exprs(msempire_eset_filtered) = dataset_msempire_matlog2 = tmp
  dataset_msempire_groupid = Biobase::pData(msempire_eset_filtered)$condition[match(colnames(dataset_msempire_matlog2), rownames(Biobase::pData(msempire_eset_filtered)))]

  msempire_example_as_tibble = msdap::import_expressionset(msempire_eset_filtered, column_fdata_protein_id = "prot.id", acquisition_mode = "dda")
  msempire_example_as_tibble$samples$condition = as.numeric(msempire_example_as_tibble$samples$group) + 1 # add metadata
  Biobase::fData(msempire_eset_filtered)$protein_id = Biobase::fData(msempire_eset_filtered)$prot.id # add metadata
  Biobase::fData(msempire_eset_filtered)$peptide_id = rownames(Biobase::fData(msempire_eset_filtered)) # add metadata
  Biobase::pData(msempire_eset_filtered)$sample_id = rownames(Biobase::pData(msempire_eset_filtered)) # add metadata


  ######
  f <- system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")
  dataset = import_dataset_skyline(f, confidence_threshold = 0.01, return_decoys = F, acquisition_mode = "dia")
  dataset = sample_metadata_custom(
    dataset,
    sample_property_regex = list(shortname = "", group = ""),
    group_regex_array = c(A = "007|009|011", B = "008|010|012")
  )
  dataset = setup_contrasts(dataset, contrast_list = list(c("A", "B")))
  dataset$samples$condition = dataset$samples$`contrast: A vs B` - 1 # add metadata

  # filter N detect in both groups
  dataset = filter_dataset(dataset,
                           filter_min_detect = 2,
                           norm_algorithm = "",
                           by_group = F,
                           all_group = T,
                           by_contrast = F)
  # optionally, select the first N peptides (that passed basic filter criteria)
  if(!is.na(subset_n_rows)) {
    valid_pepid = head(dataset$peptides %>% filter(is.finite(intensity)) %>% distinct(peptide_id) %>% pull(), subset_n_rows)
    dataset$peptides = dataset$peptides %>% filter(peptide_id %in% valid_pepid)
  }

  # convert our long-format peptide table to an ExpressionSet
  eset = tibble_as_eset(dataset$peptides %>%
                          select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, intensity=intensity_all_group) %>%
                          filter(is.finite(intensity)),
                        dataset$proteins,
                        dataset$samples)
  dataset_lfqbench_matlog2 = Biobase::exprs(eset)
  dataset_lfqbench_groupid = Biobase::pData(eset)$group[match(colnames(dataset_lfqbench_matlog2), rownames(Biobase::pData(eset)))]

  return(list(as_matrix = list(msempire_example = list(mat=dataset_msempire_matlog2, groupid=dataset_msempire_groupid),
                               lfqbench = list(mat=dataset_lfqbench_matlog2, groupid=dataset_lfqbench_groupid)),
              as_tibble = list(msempire_example = msempire_example_as_tibble,
                               lfqbench = dataset),
              as_eset = list(msempire_example = msempire_eset_filtered,
                             lfqbench = eset) ))
}

