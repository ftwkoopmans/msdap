
rm(list = ls(all.names = TRUE)) # clear everything from memory
cat("\014") # clear terminal (send the control+L character)
devtools::load_all() # load our R package

dataset = import_dataset_maxquant_evidencetxt(path = "E:/DATA/PXD007683/txt_mbr")

dataset = import_fasta(
  dataset,
  files = c(
    "E:/DATA/PXD007683/fasta/UP000002311_559292.fasta",
    "E:/DATA/PXD007683/fasta/UP000005640_9606.fasta"
  )
)

dataset = import_sample_metadata(dataset, "E:/DATA/PXD007683/oconnel_samples.xlsx")

dataset = setup_contrasts(
  dataset,
  contrast_list = list(
    c("one", "two"),
    c("one", "three"),
    c("two", "three")
  )
)

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
  dea_qvalue_threshold = 0.01,
  dea_log2foldchange_threshold = NA,
  diffdetect_min_samples_observed = 2,
  output_qc_report = TRUE,
  output_abundance_tables = TRUE,
  output_dir = "C:/temp",
  output_within_timestamped_subdirectory = TRUE,
  dump_all_data = TRUE
)

print_dataset_summary(dataset)
