
### test code; assert that the MS-DAP dataset is the same as the FragPipe summary table 'combined_peptide.tsv'   (ideally we'd just parse that, but it lacks retention-time data)

#
# dataset = import_dataset_fragpipe_ionquant(path = "C:/DATA/PXD007683_oconnel-2018/fragpipe_output", acquisition_mode = "dda", collapse_peptide_by = "")
#
# ref = readr::read_tsv("C:/DATA/PXD007683_oconnel-2018/fragpipe_output/combined_peptide.tsv")
#
# sample_idmap = c(
#   "a05191"="one_1",
#   "a05192"="one_2",
#   "a05194"="one_3",
#   "a05195"="two_1",
#   "a05196"="two_2",
#   "a05197"="two_3",
#   "a05198"="two_4",
#   "a05199"="three_1",
#   "a05200"="three_2",
#   "a05201"="three_3",
#   "a05202"="three_4"
# )
#
# tib_plot_msdap = dataset$peptides %>%
#   group_by(sample_id, sequence_plain) %>%
#   summarise(
#     protein_id = protein_id[1],
#     intensity = log2(sum(2^intensity)),
#     .groups = "drop"
#   ) %>%
#   left_join(tibble::enframe(sample_idmap, name = "sample_id", value = "sample_label"), by = "sample_id")
#
# tib_plot_ref = ref %>%
#   select(sequence_plain = `Peptide Sequence`, protein_id_ref = `Protein ID`, !!grep(" intensity$", colnames(ref), value=T, ignore.case = T)) %>%
#   pivot_longer(cols = -c("sequence_plain", "protein_id_ref"), names_to = "sample_id", values_to = "intensity") %>%
#   filter(is.finite(intensity) & intensity > 0) %>%
#   mutate(
#     intensity = log2(intensity),
#     sample_id = gsub(" intensity$", "", sample_id, ignore.case = T)
#   )
#
# # merged table with peptide-level data from both data sources combined
# tib_plot = tib_plot_msdap %>%
#   left_join(tib_plot_ref %>% rename(intensity_ref = intensity), by = c("sample_label"="sample_id", "sequence_plain"="sequence_plain"))
#
# # check 1: peptide abundances are correlated (some variance might apply, e.g. due to normalization)
# print( ggplot(tib_plot, aes(intensity_ref, intensity)) +
#          geom_point() +
#          geom_abline(intercept = 0, slope = 1, col = "red") +
#          facet_wrap(.~sample_id) )
#
# # check 2: protein identifier is the same for all peptides
# print( table(gsub(";.*", "", tib_plot$protein_id) == tib_plot$protein_id_ref) )
#
#
#
#
#
#
# # expect error; parameter errors
# dataset = import_dataset_fragpipe_ionquant("D:/DOESNOTEXIST", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "test", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "test")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified", confidence_threshold = -1)
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified", confidence_threshold = "0")
#
# # expect error; no metadata thus import error
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# # expect warning; metadata is not unique per sample thus input is assumed to be fractionated (multiple sample/rawfile per 'experiment')
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experiment1value-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# unique(dataset$peptides$sample_id) # show combined sample_id
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# # expect success; valid datasets
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentunique-nobiorep-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentunique-nobiorep-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-nombr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_modified")
#
# # expect warning; collapse param = keep modseq*charge
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "")
# # expect succes; now merging by plain sequence
# dataset = import_dataset_fragpipe_ionquant("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_experimentproper-biorepproper-mbr", acquisition_mode = "dda", collapse_peptide_by = "sequence_plain")
#
#
#
#
# dataset = import_dataset_fragpipe_psm_file("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata/psm.tsv", acquisition_mode = "dda", collapse_peptide_by = "")
# dataset$peptides
# dataset1 = import_dataset_fragpipe_psm_file("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata/psm.tsv", intensity_sum = T, acquisition_mode = "dda", collapse_peptide_by = "")
# dataset2 = import_dataset_fragpipe_psm_file("D:/DATA/fragpipe_test_oconnel2018/fragpipe_18_no-experiment-metadata/psm.tsv", intensity_sum = F, acquisition_mode = "dda", collapse_peptide_by = "")
#
# dataset1$peptides
# dataset2$peptides
# hist(dataset1$peptides$intensity - dataset2$peptides$intensity)
# table(dataset1$peptides$intensity != dataset2$peptides$intensity)



