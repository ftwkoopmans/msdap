
rm(list = ls(all.names = TRUE)) # clear everything from memory
cat("\014") # clear terminal (send the control+L character)
devtools::load_all() # load our R package

dataset = import_dataset_metamorpheus(path = "C:/VU/code/R/msdap/docker/temp/exampledata/dataset_Klaassen2018_pmid26931375", protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified")

dataset = import_fasta(
  dataset,
  files = c("C:/VU/fasta/UniProt_2018-05/UP000000589_10090.fasta",
            "C:/VU/fasta/UniProt_2018-05/UP000000589_10090_additional.fasta")
)

dataset = remove_proteins_by_name(dataset, regular_expression = "ig \\S+ chain|keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)") # optional

dataset = sample_metadata_custom(
  dataset,
  # a list of regular expressions to extract a short name and a group label from each full sample name
  sample_property_regex = list(
    shortname = c(".*CB_([A-Z]+.\\d).*", "\\1"), # regex dictates what part of the filename should be removed
    group = c(".*CB_([A-Z]+)_.*", "\\1"),
    fraction = c(".*\\D(\\d+)\\Dqtof.*", "\\1")
  ),
  # optionally, sort the groups in specified order
  group_order = c("WT", "KO")
)
dataset = setup_contrasts(dataset, contrast_list = list(c("WT", "KO")))

dataset = analysis_quickstart(
  dataset,
  # filter_min_detect = 1, # !! many proteins are only found with MBR in the knockout
  filter_min_detect = 0,
  filter_min_quant = 3,
  filter_fraction_detect = 0,
  filter_fraction_quant = 0.75,
  filter_by_contrast = TRUE,
  filter_topn_peptides = 0,
  filter_min_peptide_per_prot = 1,
  norm_algorithm = c("vsn", "modebetween_protein"),
  dea_algorithm = c("deqms", "msempire", "msqrob"),
  dea_qvalue_threshold = 0.05,
  dea_log2foldchange_threshold = NA,
  diffdetect_min_samples_observed = 2,
  output_qc_report = TRUE,
  output_abundance_tables = TRUE,
  output_dir = "C:/temp",
  output_within_timestamped_subdirectory = TRUE,
  dump_all_data = TRUE
)

print_dataset_summary(dataset)
