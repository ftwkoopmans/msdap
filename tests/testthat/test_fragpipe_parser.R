#
#
# ### test: import various "datasets" generated from the same raw data, by processing them with various FragPipe versions and settings in the workflow tab ("experimental designs")
#
# devtools::load_all()
#
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_15.0_experiment=filename_biorep=empty_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_15.0_experiment=sample_biorep=empty_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_15.0_experiment=condition_biorep=replicate_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_18.0_experiment=filename_biorep=empty_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_18.0_experiment=sample_biorep=empty_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_18.0_experiment=condition_biorep=replicate_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_20.0_experiment=filename_no-biorep_msfragger-default-settings_extra-mod-phospho", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_22.0_experiment=filename_biorep=empty_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_22.0_experiment=sample_biorep=empty_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant(path = "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/fragpipe/version_22.0_experiment=condition_biorep=replicate_collapse-fractions_msfragger-default-settings", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
#
#
# ### now apply the full MS-DAP pipeline and verify that we obtain expected results for this study with mostly known true-postives
#
# dataset = import_fasta(dataset, "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/uniprot_2024_03/2024-07-28-decoys-contam-uniprot_2024_03_UP000000589_10090_full_proteome.fasta.fas")
# dataset = remove_proteins_by_name(dataset, regular_expression = "ig \\S+ chain|keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)") # optional
#
# dataset = import_sample_metadata(dataset, "C:/VU/projects/Remco - Shisa6 IP WT-KO (2015)/mysamples.xlsx")
# dataset = setup_contrasts(dataset, contrast_list = list(c("wt", "ko")))
#
# dataset = analysis_quickstart(
#   dataset,
#   # filter_min_detect = 1, # !! many proteins are only found with MBR in the knockout
#   filter_min_detect = 0,
#   filter_min_quant = 3,
#   filter_fraction_detect = 0,
#   filter_fraction_quant = 0.75,
#   filter_by_contrast = FALSE,
#   filter_topn_peptides = 0,
#   filter_min_peptide_per_prot = 1,
#   norm_algorithm = c("vsn", "modebetween_protein"),
#   dea_algorithm = c("deqms", "msempire", "msqrob"),
#   dea_qvalue_threshold = 0.05,
#   dea_log2foldchange_threshold = NA,
#   diffdetect_min_samples_observed = 2,
#   output_qc_report = TRUE,
#   output_abundance_tables = TRUE,
#   output_dir = "C:/temp",
#   output_within_timestamped_subdirectory = FALSE,
#   dump_all_data = TRUE
# )
#
# print_dataset_summary(dataset)
#
#
# ###
# hgnc_table = hgnc_lookuptable("C:/DATA/HGNC/hgnc_complete_set__2024-01-01.txt")  # <<EDIT THIS FILENAME>>
# # MGI Mouse database, only needed for Mouse datasets
# mgi_table = mgi_lookuptable("C:/DATA/mgi/MRK_SwissProt_TrEMBL__2024-07-13.rpt")  # <<EDIT THIS FILENAME>>
# # RGD Rat database, only needed for Rat datasets
# # rgd_table = rgd_lookuptable("C:/DATA/rgd/GENES_RAT__2024-07-13.txt")  # <<EDIT THIS FILENAME>>
#
# plot_differential_detect(dataset, zscore_threshold = 7)
#
# # refer to this function's R documentation for additional details
# tmp = export_stats_genesummary(
#   dataset,
#   # set NA to ignore differential detection (default), or define an absolute zscore threshold (e.g. 6)
#   # ! when using diffdetect, first review the zscore distributions using MS-DAP function: plot_differential_detect()
#   diffdetect_zscore_threshold = 6,
#   # options for dealing with proteingroups that map to multiple genes:
#   # "leading_gene" = use only the leading/first gene per proteingroup. e.g. "GRIA1;GRIA2" is mapped to "GRIA1" and will overwrite a specific "GRIA1" proteingroup if (and only if) the former has a lower p-value
#   # "prio_specific" = discard ambiguous proteingroups if their leading/first gene overlaps with a specific proteingroup. e.g. "GRIA1;GRIA2" is mapped to "GRIA1" only if there is no specific/unambiguous "GRIA1" proteingroup
#   # "only_specific" = all ambiguous proteingroups are discarded
#   gene_ambiguity = "prio_specific", # prio_specific and only_specific are recommended
#   # HGNC data table is always required
#   hgnc = hgnc_table,
#   # example: for Mouse datasets, provide the MGI database for more accurate ID mapping
#   xref = mgi_table,
#   # example: for Rat datasets, provide the RGD database for more accurate ID mapping
#   # xref = rgd_table,
#   # write output files to the same directory as main MS-DAP results by using the `dataset$output_dir` variable
#   # alternatively, provide a valid path. e.g. output_dir="C:/temp" # use forward slashes
#   output_dir = "C:/temp" #dataset$output_dir
# )
#
# tmp |> filter(dea_algorithm == "msqrob") |> print(n=50)
# # dataset$de_proteins |> left_join(dataset$proteins |> select(protein_id, gene_symbols_or_id)) |> filter(grepl("GRIA1", gene_symbols_or_id) & dea_algorithm == "deqms")
# # tmp |> filter(grepl("GRIA|CACNG|CNIH|PRRT", gene_symbols_or_id) & dea_algorithm == "msqrob")
# # dataset$peptides |> left_join(dataset$proteins |> select(protein_id, gene_symbols_or_id)) |> filter(grepl("CACNG", gene_symbols_or_id)) |> select(gene_symbols_or_id, sample_id, peptide_id, intensity, detect) |> arrange(gene_symbols_or_id, sample_id, peptide_id)
# # dataset$dd_proteins |> left_join(dataset$proteins |> select(protein_id, gene_symbols_or_id)) |> filter(grepl("GRIA|CACNG", gene_symbols_or_id) & type == "detect") |> arrange(gene_symbols_or_id, desc(abs(zscore)))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ### another set of FragPipe searches for testing the import function
#
# # expect errors due to invalid parameters
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("D:/DOESNOTEXIST", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified") )
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified") )
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "test", collapse_peptide_by = "sequence_modified") )
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "test") )
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified", confidence_threshold = -1) )
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified", confidence_threshold = "0") )
#
# # expect errors @ edge-case; no metadata thus import error
# testthat::expect_error( dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified") )
#
# # edge-case; metadata is not unique per sample thus input is assumed to be fractionated (multiple sample/rawfile per 'experiment')
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# unique(dataset$peptides$sample_id) # show combined sample_id
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# # typical FragPipe configurations
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentunique-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentunique-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
# # expect warning; collapse param = keep modseq*charge
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "")
# # expect succes; now merging by plain sequence
# dataset = import_dataset_fragpipe_ionquant("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_plain")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ### expect warning; import PSM file directly
# dataset1 = import_dataset_fragpipe_psm_file("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata/psm.tsv", intensity_sum = T, acquisition_mode = "dda", collapse_peptide_by = "")
# dataset2 = import_dataset_fragpipe_psm_file("E:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata/psm.tsv", intensity_sum = F, acquisition_mode = "dda", collapse_peptide_by = "")
#
# # now plot the "highest intensity PSM" approach versus the "sum of all PSM intensities"
# tmp = bind_rows(dataset1$peptides |> select(sample_id, peptide_id, intensity) |> mutate(id = paste(sample_id, peptide_id), config="sum"),
#                 dataset2$peptides |> select(sample_id, peptide_id, intensity) |> mutate(id = paste(sample_id, peptide_id), config="top1")) |>
#   pivot_wider(id_cols = c("sample_id", "id"), names_from = "config", values_from = "intensity")
#
# ggplot(tmp, aes(sum, top1)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, colour = "red") +
#   facet_wrap(.~sample_id)
