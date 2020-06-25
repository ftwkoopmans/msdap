rm(list = ls(all.names = TRUE)) # clear everything from memory
cat("\014") # clear terminal (send the control+L character)
devtools::load_all() # load our R package


filter_topn_peptides = 5

list_datasets = list()

############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################

ds_lfqbench = import_dataset_spectronaut(
  # file_spectronaut = "C:/Users/Frank/Downloads/Spectronaut_20150522_163227_TTOF5600_64w_newLib_Report.tsv", # original report
  filename = "C:/DATA/lfqbench/20191017_130053_lfqbench_5600_64var_Report.xls", # our re-run in spectronaut 13
  confidence_threshold = 0.01,
  use_normalized_intensities = F,
  do_plot = F
)
ds_lfqbench = setup_contrasts(sample_metadata_custom(ds_lfqbench, sample_property_regex = list(shortname = "", group = ""), group_regex_array = c(A = "007|009|011", B = "008|010|012")) , contrast_list = list(c("A", "B")))

list_datasets[[length(list_datasets) + 1]] = list(dataset = ds_lfqbench,
                                                  name = "lfqbench_spectronaut13-nonorm_mindetect2",
                                                  filter_min_detect = 2,
                                                  filter_fraction_detect = 0,
                                                  filter_min_quant = 3,
                                                  filter_fraction_quant = 0,
                                                  filter_min_peptide_per_prot = 1,
                                                  filter_topn_peptides = filter_topn_peptides,
                                                  protein_classification_regex = c(background="_HUMA", foreground="_YEAS", discard="_ECOL"))
list_datasets[[length(list_datasets) + 1]] = list(dataset = ds_lfqbench,
                                                  name = "lfqbench_spectronaut13-nonorm_mindetect3",
                                                  filter_min_detect = 3,
                                                  filter_fraction_detect = 0,
                                                  filter_min_quant = 3,
                                                  filter_fraction_quant = 0,
                                                  filter_min_peptide_per_prot = 1,
                                                  filter_topn_peptides = filter_topn_peptides,
                                                  protein_classification_regex = c(background="_HUMA", foreground="_YEAS", discard="_ECOL"))




############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################

ds_oconnel = import_dataset_metamorpheus(path = "C:/DATA/PXD007683/2019-10-14-19-20-46/Task2-SearchTask", protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified")
ds_oconnel = import_fasta(ds_oconnel, files = c("C:/DATA/PXD007683/fasta/UP000002311_559292.fasta", "C:/DATA/PXD007683/fasta/UP000005640_9606.fasta"))
ds_oconnel = setup_contrasts(sample_metadata_custom(ds_oconnel,group_regex_array = c(one="a05191|a05192|a05194", two="a05195|a05196|a05197|a05198", three="a05199|a05200|a05201|a05202"),
                                                    sample_property_regex = list(shortname = "_.*", group = ""),
                                                    group_order = c("one", "two", "three")),
                             contrast_list = list(c("one", "three"), c("one", "two"), c("two", "three")))


list_datasets[[length(list_datasets) + 1]] = list(dataset = ds_oconnel,
                                                  name = "OConnel2018_metamorpheus",
                                                  filter_min_detect = 1,
                                                  filter_fraction_detect = 0,
                                                  filter_min_quant = 3,
                                                  filter_fraction_quant = 0,
                                                  filter_min_peptide_per_prot = 1,
                                                  filter_topn_peptides = filter_topn_peptides,
                                                  protein_classification_regex = c(background="OX=9606", foreground="OX=559292", discard=""))




############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################

ds_lim = import_dataset_maxquant_evidencetxt("C:/DATA/PXD014415/txt", collapse_peptide_by = "sequence_modified", remove_shared = T)
ds_lim = sample_metadata_custom(ds_lim, sample_property_regex = list(
  shortname = "_90min.*", # regex dictates what part of the filename should be removed
  group = "(^[[:alnum:]]+_)|(_90min.*)"))



### manipulate samples; only keep the human samples, then split into two groups.
# Sample numbers/IDs might be the measurement order, so we assign alternating groups to negate potential batch effects
ds_lim$samples = ds_lim$samples %>% filter(group == "human")
ds_lim$samples$group = LETTERS[1:2][1:nrow(ds_lim$samples) %% 2 + 1]
# establish a contrast for downstream statistics; compare group A vs group B
ds_lim = setup_contrasts(ds_lim, contrast_list = list(c("A", "B")))

# first 8 samples in each group
samples_id_A = ds_lim$samples %>% filter(group == "A") %>% pull(sample_id) %>% head(n=8)
samples_id_B = ds_lim$samples %>% filter(group == "B") %>% pull(sample_id) %>% head(n=8)

# remove peptide detected in less than 2 replicates within either group
pep_mindetect_2 = inner_join(ds_lim$peptides %>% filter(sample_id %in% samples_id_A & detect==T) %>% count(peptide_id, name = "n_grpA") %>% filter(n_grpA>=2),
                             ds_lim$peptides %>% filter(sample_id %in% samples_id_B & detect==T) %>% count(peptide_id, name = "n_grpB") %>% filter(n_grpB>=2),
                             by="peptide_id")

# apply filter to the peptide data table
ds_lim$peptides = ds_lim$peptides %>% filter(sample_id %in% c(samples_id_A, samples_id_B) & peptide_id %in% pep_mindetect_2$peptide_id)



### manipulate peptide abundances;
# note: some will be impossible to detect downstream because there are not enough data points in either group
# will also happen in real datasets. Only working on peptides/proteins with many data points is not representative of a typical workflow
mock_foldchange = 1.2

# flag 10% of unique proteins at random to set artificial fold-change (between the artifical groups)
uprot = ds_lim$peptides %>% distinct(protein_id) %>% pull()
uprot_diff = uprot[1:floor(0.1 * length(uprot))] # as a random set, just take first N protein IDs
ds_lim$peptides$is_diff = ds_lim$peptides$protein_id %in% uprot_diff
# update protein descriptions so we can recognize them as 'spike in, should be different' downstream
ds_lim$proteins = ds_lim$proteins %>% filter(protein_id %in% uprot) %>% mutate(is_diff = protein_id %in% uprot_diff)
ds_lim$proteins$fasta_headers[ds_lim$proteins$is_diff] = paste0("DIFF_", ds_lim$proteins$fasta_headers[ds_lim$proteins$is_diff])
# debug; table(ds_lim$proteins$is_diff)

# increase 'spike-in' peptide abundances by x%
rows = ds_lim$peptides$is_diff & ds_lim$peptides$sample_id %in% samples_id_B
ds_lim$peptides$intensity[rows] = ds_lim$peptides$intensity[rows] + log2(mock_foldchange)


####


for(n_replicates in c(8,6,4)) { #n_replicates=4
  valid_sample_ids = c(ds_lim$samples %>% filter(group == "A") %>% pull(sample_id) %>% head(n=n_replicates),
                       ds_lim$samples %>% filter(group == "B") %>% pull(sample_id) %>% head(n=n_replicates))

  ds = ds_lim
  ds$peptides = ds$peptides %>% filter(sample_id %in% valid_sample_ids)
  ds$samples = ds$samples %>% filter(sample_id %in% valid_sample_ids)

  list_datasets[[length(list_datasets)+1]] = list(dataset = ds,
                                                  name = sprintf("DDA_in-silico_spike-in_nrep-%s_FC-%s", n_replicates, mock_foldchange),
                                                  filter_min_detect = 1,
                                                  filter_fraction_detect = .25,
                                                  filter_min_quant = 3,
                                                  filter_fraction_quant = .75,
                                                  filter_min_peptide_per_prot = 1,
                                                  filter_topn_peptides = filter_topn_peptides,
                                                  protein_classification_regex = c(background="^(?!DIFF)", foreground="^DIFF", discard="")) # negative lookahead regex
}



############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
# save(list_datasets, file = "C:/temp/benchmark_datasets.RData", compress = T)



cl <<- initialize_multiprocessing()

####
roc_data = NULL

for(ds in list_datasets) { # ds=list_datasets[[4]]
  # classify proteins
  ds$dataset$proteins$classification = regex_classification(ds$dataset$proteins$fasta_headers, regex=ds$protein_classification_regex)
  print(table(ds$dataset$proteins$classification))


  # #######################################################################################################################################################################
  var_overall = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
    normalize_vwmb(x_as_log2, groups=NA, metric_within="var")
  }
  # N_var_a = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=NA, metric_within="var", exclude_rare_detections = F, exclude_most_variation = F)
  # }
  # N_var_a_exdc = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=NA, metric_within="var", exclude_rare_detections = T, exclude_most_variation = T)
  # }
  # N_var_mode = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=mask_sample_groups, metric_within="var", metric_between = "mode", exclude_rare_detections = F, exclude_most_variation = F)
  # }
  # N_var_mode_exd = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=mask_sample_groups, metric_within="var", metric_between = "mode", exclude_rare_detections = T, exclude_most_variation = F)
  # }
  # N_var_mode_exdc = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=mask_sample_groups, metric_within="var", metric_between = "mode", exclude_rare_detections = T, exclude_most_variation = T)
  # }
  #
  # N_modebetween = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=mask_sample_groups, metric_within="", metric_between = "mode", exclude_rare_detections = F, exclude_most_variation = F)
  # }
  # N_modebetween_exdc = function(x_as_log2, mask_sample_groups = NA, x_mask_values_to_use = NA) {
  #   mynorm(x_as_log2, groups=mask_sample_groups, metric_within="", metric_between = "mode", exclude_rare_detections = T, exclude_most_variation = T)
  # }
  # #######################################################################################################################################################################


  ##################### prep data #####################
  # for(norm_algorithm in list("vwmb")) { # norm_algorithm = "vwmb"
  for(norm_algorithm in list("vwmb", "msempire", "var_overall", "loess", "vsn", "rlr", c("rlr","modebetween"), c("vsn","modebetween"), c("loess","modebetween"))) { # norm_algorithm = "vwmb"
    # for(norm_algorithm in list("mode", "msempire", "N_var_a", "N_var_a_exdc", "N_var_mode", "N_var_mode_exd", "N_var_mode_exdc", "loess", "vsn", c("vsn","N_modebetween"), c("vsn","N_modebetween_exdc"))) { # norm_algorithm = "N_var_test2"
    # for(norm_algorithm in c("mode", "msempire", "vsn", "loess", "N_var_a", "N_var_wb", "N_var_mode")) { # norm_algorithm = "mode"
    # for(norm_algorithm in c("mode", "msempire", "N_var_a", "N_var_wb", "N_var_mode", "N_var_test1", "N_var_test2", "N_var_test3", "N_var_test4")) { # norm_algorithm = "N_var_test2"
    # for(norm_algorithm in c("mode", "msempire", "N_var_wb", "N_var_mode", "N_var_wb_no", "N_var_mode_no")) { # norm_algorithm = ""

    cat("***", ds$name, norm_algorithm, "***\n")

    # filter
    ds$dataset = filter_dataset(ds$dataset,
                                filter_min_detect = ds$filter_min_detect, filter_min_quant = ds$filter_min_quant,
                                filter_fraction_detect = ds$filter_fraction_detect, filter_fraction_quant = ds$filter_fraction_quant,
                                filter_min_peptide_per_prot = ds$filter_min_peptide_per_prot, filter_topn_peptides = ds$filter_topn_peptides,
                                norm_algorithm = norm_algorithm,
                                by_group = F, all_group = F, by_contrast = T) # must do by contrast, because msempire normalization does not support more than 2 groups (eg; oconnel dataset with 3 groups)
    # dea
    ds$dataset = dea(ds$dataset, qval_signif = 0.01, algo_de = c("ebayes", "msempire", "msqrob"))
    # ds$dataset = dea(ds$dataset, qval_signif = 0.01, algo_de = c("ebayes", "msempire"))
    # ds$dataset = dea(ds$dataset, qval_signif = 0.01, algo_de = "ebayes")
    ## print number signif: ds$dataset$de_proteins %>% dplyr::filter(signif) %>% dplyr::count(algo_de, contrast); dplyr::n_distinct(ds$dataset$de_proteins$protein_id)
    ## debug normalization:
    # x = as_matrix_except_first_column(ds$dataset$peptides %>% pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "intensity_all_group"))
    # boxplot(x,las=2); y = N_var_test4(x, mask_sample_groups = ds$dataset$samples$group[match(colnames(x), ds$dataset$samples$sample_id)]); boxplot(y,las=2)

    # add protein metadata  +  add minlog10 qvalue
    ds$dataset$de_proteins = ds$dataset$de_proteins %>%
      left_join(ds$dataset$proteins, by="protein_id") %>%
      mutate(predictor = classification) %>%
      group_by(algo_de, contrast) %>%
      mutate(qvalue_minlog10 = minlog10(qvalue)) %>%
      ungroup()

    # store ROC data; only foreground/background proteins  +  add this loop's normalization
    roc_data = bind_rows(roc_data, ds$dataset$de_proteins %>%
                           filter(predictor %in% c("background", "foreground")) %>%
                           add_column(algo_norm = paste(norm_algorithm, collapse="&"), ds_name=ds$name) )
  }
}

parallel::stopCluster(cl)


save(list_datasets, roc_data, file = sprintf("C:/temp/benchmark_datasets_topN-%s.RData", filter_topn_peptides), compress = T)

