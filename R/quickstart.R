#' MS-DAP: A package for Downstream Analysis of Proteomics.
#'
#' http://github.com/ftwkoopmans/msdap
#'
#' @docType package
#' @name msdap
#'
# some additional imports to get rid of devtools::check() warnings
#' @import ggplot2
#' @importFrom grDevices colorRampPalette dev.off graphics.off pdf rainbow
#' @importFrom graphics abline boxplot legend lines mtext par plot points text
#' @importFrom stats .lm.fit contrasts cov density ecdf mad median na.exclude na.omit p.adjust p.adjust.methods pt quantile resid residuals rt sd sigma weights loess optim predict
#' @importFrom utils capture.output combn data head relist tail help packageVersion
#' @importFrom rlang .data :=
NULL



#' Quickstart for analyses in this pipeline
#'
#' all-in-one function that covers the vast majority of use-cases of analyzing a dataset imported into MS-DAP.
#' (assuming you already loaded peptide data, sample metadata and fasta files using MS-DAP import functions).
#'
#' @section Filtering:
#'
#' Peptide filter criteria applied to replicate samples within a sample group.
#' params; filter_min_detect, filter_fraction_detect, filter_min_quant, filter_fraction_quant.
#' You only have to provide active filters (but specify at least 1), filters/settings you do not specify don't do anything by default.
#'
#' Settings:
#' for DDA: at least 1~2 detect (MS/MS ID) and quantified in at least ~75% of replicates.
#' for DIA: detect (confidence score < threshold) in at least ~75% of replicates (because for DIA, you typically have an abundance value in each sample regardless of the identifier confidence score).
#' If there are only 3 replicates, we recommend filtering such that there are at least 3 datapoints to work with in differential expression analysis.
#'
#' Taken together, recommended settings for a DDA dataset with 3~8 replicates in each sample group look like this;
#'
#' \code{
#' filter_min_detect = 1 (or zero to fully rely on MBR), filter_fraction_detect = 0.25 (or zero to fully rely on MBR), filter_min_quant = 3, filter_fraction_quant = 0.75
#' }
#'
#' Analogous for DIA;
#'
#' \code{
#' filter_min_detect = 3, filter_fraction_detect = 0.75
#' }
#'
#' @section Filter within contrast vs using all groups:
#'
#' Two distinct approaches to selecting peptides can be used for differential expression analysis: 1) 'within contrast' and 2) 'apply filter to all sample groups'.
#'
#' 1) Determine within each contrast (eg; group A vs group B) what peptides can be used by applying above peptide filter criteria and then apply normalization to this data subset.
#' Advantaguous in datasets with many groups; this maximizes the number of peptides used in each contrast (eg; let peptide p be observed in groups A and B, not in C. we'd want to use it in A vs B, not in A vs C).
#' As a disadvantage, this complicates interpretation since the exact data used is different in each contrast (slightly different peptides and normalization in each contrast).
#'
#' 2) Apply above filter criteria to each sample group (eg; a peptide must past these filter rules in every sample group) and then apply normalization
#'
#' This data matrix is then used for all downstream statistics
#'
#' Advantage; simple and robust
#'
#' Disadvantage; potentially miss out on (group-specific) peptides/data-points that may fail filter criteria in just 1 group, particularly in large datasets with 4+ groups
#'
#' Set \code{filter_within_contrast = FALSE} for this option
#'
#' Note; if there are just 2 sample groups (eg; WT vs KO), this point is moot as both approaches are the same
#'
#' @section Normalization:
#' normalization algorithms are applied to the peptide-level data matrix.
#' options: "" (empty string disables normalization), "vsn", "loess", "rlr", "msempire", "vwmb", "modebetween", "modebetween_protein" (this balances foldchanged between sample groups. Highly recommended, see MS-DAP manuscript)
#'
#' You can combine normalizations by providing an array of options to apply subsequential normalizations.
#'
#' For instance, \code{norm_algorithm = c("vsn", "modebetween_protein")} applies the vsn algorithm (quite strong normalization reducing variation) and then balances between-group protein-level foldchanges with modebetween normalization.
#'
#' Benchmarks have shown that c("vwmb", "modebetween_protein") and c("vsn", "modebetween_protein") are the optimal strategies, see MS-DAP manuscript.
#'
#' @section Differential Expression Analysis:
#'
#' Statistical models for differential expression analysis
#'
#' MSqRob is recommended for most cases; a peptide-level model that is highly sensitive and quite robust. Reference: https://github.com/statOmics/MSqRob
#'
#' MS-EmpiRe a peptide-level model that works especially well for DDA data. Reference: https://github.com/zimmerlab/MS-EmpiRe
#'
#' eBayes is robust but conservative, using the limma package to apply moderated t-tests on protein-level abundances. Reference: https://doi.org/doi:10.18129/B9.bioc.limma
#'
#' options: ebayes, deqms, msempire, msqrob, msqrobsum
#'
#' You can simply apply multiple DEA models in parallel by supplying an array of options. The output of each model will be visualized in the PDF report and data included in the output Excel report.
#' e.g.; \code{dea_algorithm = c("ebayes", "deqms", "msempire", "msqrob")}
#'
#'
#'
#' @param dataset a valid dataset object generated upstream by an MS-DAP import function. For instance, import_dataset_skyline() or import_dataset_maxquant_evidencetxt()
#' @param filter_min_detect in order for a peptide to 'pass' in a sample group, in how many replicates must it be detected?
#' @param filter_fraction_detect in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be detected?
#' @param filter_min_quant in order for a peptide to 'pass' in a sample group, in how many replicates must it be quantified?
#' @param filter_fraction_quant in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be quantified?
#' @param filter_min_peptide_per_prot in order for a peptide to 'pass' in a sample group, how many peptides should be available after detect filters? 1 is default, but 2 can be a good choice situationally (eg; to not rely on proteins with just 1 quantified peptide)
#' @param filter_topn_peptides maximum number of peptides to maintain for each protein (from the subset that passes above filters, peptides are ranked by the number of samples where detected and their variation between replicates).
#' @param filter_by_contrast should the above filters be applied to all sample groups, or only those tested within each contrast? Enabling this optimizes available data in each contrast, but increases the complexity somewhat as different subsets of peptides are used in each contrast and normalization is applied separately.
#' @param norm_algorithm normalization algorithms. Can be empty string to skip normalization, or any algorithm listed by msdap::normalization_algorithms(). Provide an array of options to run each algorithm consecutively, for instance; c("vsn", "modebetween_protein") to first apply vsn normalization and then correct between-group ratios such that the protein-level log2-foldchange mode is zero
#' @param dea_protein_rollup strategy for combining peptides to proteins as used in DEA algorithms that first combine peptides to proteins and then apply statistics, like ebayes and deqms. options: maxlfq, sum. The former applies the MaxLFQ algorithm, the latter employs the 'classic' strategy of summing all peptides per protein. See further rollup_pep2prot()
#' @param dea_algorithm algorithms for differential expression analysis. options: ebayes, deqms, msqrobsum, msempire, msqrob (to run multiple, provide an array)
#' @param dea_qvalue_threshold threshold for significance of adjusted p-values in figures and output tables. Output tables will also include all q-values as-is
#' @param dea_log2foldchange_threshold threshold for significance of log2 foldchanges. Set to zero to disregard or a positive value to apply a cutoff to absolute log2 foldchanges. MS-DAP can also perform a bootstrap analyses to infer a reasonable threshold by setting this parameter to NA
#' @param diffdetect_min_samples_observed for differential detection only; minimum number of samples where a protein should be observed at least once by any of its peptides (in either group) when comparing a contrast of group A vs B
#' @param pca_sample_labels whether to use sample names or a numeric ID as labels in the PCA plot. options: "auto" (let the code decide, default), "shortname" (use sample shortnames), "index" (auto-generated numeric ID), "index_asis" (same as index option and specifically disable label overlap reduction)
#' @param output_dir output directory where all output files should be stored. If not an existing directory, it will be created
#' @param output_within_timestamped_subdirectory optionally, automatically create a subdirectory (within output_dir) that has the current date&time as name and store results there. options: FALSE, TRUE
#' @param output_abundance_tables whether to write protein-level data matrices to file. options: FALSE, TRUE
#' @param output_qc_report whether to create the Quality Control report. options: FALSE, TRUE . Highly recommended to set to TRUE (default). Set to FALSE to skip the report PDF (eg; to only do differential expression analysis and skip the time-consuming report creation)
#' @param dump_all_data if you're interested in performing custom bioinformatic analyses and want to use any of the data generated by this tool, you can dump all intermediate files to disk. Has performance impact so don't enable by default. options: FALSE, TRUE
#' @importFrom parallel stopCluster
#' @importFrom openxlsx write.xlsx
#' @importFrom data.table fwrite
#' @export
analysis_quickstart = function(dataset,
                               # peptide filter criteria applied within each sample group
                               filter_min_detect = 0, filter_fraction_detect = 0, filter_min_quant = 0, filter_fraction_quant = 0,
                               # respective criteria on protein level
                               filter_min_peptide_per_prot = 1, filter_topn_peptides = 0,
                               # apply filter to each sample group, or only apply filter within relevant sample groups being compared in a contrast?
                               filter_by_contrast = FALSE,
                               # normalization algorithms for peptide-level data. Options listed by msdap::normalization_algorithms() , provide an array of options to run multiple sequentially
                               norm_algorithm = c("vwmb", "modebetween_protein"),
                               # DEA algorithms
                               dea_protein_rollup = "maxlfq",
                               dea_algorithm = c("deqms", "msqrob"), # options; ebayes, deqms, msempire, msqrob, msqrobsum (to run multiple in parallel, provide an array of options)
                               # define thresholds for significant proteins
                               dea_qvalue_threshold = 0.01,
                               dea_log2foldchange_threshold = 0, # if NA, infer from bootstrap
                               # differential detection
                               diffdetect_min_samples_observed = 3,
                               # plot options
                               pca_sample_labels = "auto",
                               # output data
                               output_abundance_tables = TRUE,
                               output_qc_report = TRUE,
                               output_dir,
                               output_within_timestamped_subdirectory = TRUE,
                               dump_all_data = FALSE) {

  # convert to absolute paths. example; user used setwd() and now sets "output" as the output_dir
  # only forward slashes and remove redundant. normalizePath() also cleans slashes etc. like our path_clean_slashes (which does more)
  output_dir = normalizePath(paste0(output_dir, "/"), winslash = "/", mustWork = FALSE)
  # create output directory if it does not exist
  if(!dir.exists(output_dir)) {
    append_log(paste("output directory does not exist yet, creating;", output_dir), type = "info")
    dir.create(output_dir, recursive = T)
    if(!dir.exists(output_dir)) {
      append_log(paste("failed to create output directory;", output_dir), type = "error")
    }
  }

  if(file.access(output_dir, mode = 2) != 0) {
    append_log(paste("no write access to the output directory;", output_dir), type = "error")
  }

  if(filter_min_detect == 0 && filter_fraction_detect == 0 && filter_min_quant == 0 && filter_fraction_quant == 0) {
    append_log("must specify at least one parameter for filtering peptides. Any of; filter_min_detect, filter_fraction_detect, filter_min_quant, filter_fraction_quant", type = "error")
  }

  ##### prior to any analyses, check if the data is fractionated. If so, merge all sample fractions prior to downstream analyses
  if("fraction" %in% colnames(dataset$samples)) {
    append_log("sample metadata contains samples with multiple fractions, 'sample_id' with the same 'shortname' are now merged by summation of their respective peptide intensities", type = "info")
    dataset = merge_fractionated_samples(dataset)
  }


  ### check input
  check_dataset_integrity(dataset)
  # extra check; don't include any decoys in downstream data analysis
  if(any(dataset$peptides$isdecoy)) {
    append_log("peptides tibble must contain no decoy entries at this point (all values in 'isdecoy' column must be FALSE). For typical workflows, make sure to set return_decoys=FALSE when importing data.", type = "error")
  }


  ### prior to starting computation, check if standard output files are available
  if(output_within_timestamped_subdirectory) {
    # If there is an environment variable that hardcodes the timezone (an integer offset from UTC), use it to adjust the timezone
    # We use it to control timezone differences between the Docker host and the Docker container
    # (windows host systems uses different timezone abbreviations than unix. This solution is straight forward to implement and does introduce additional dependencies)
    # Windows PowerShell: [System.TimeZone]::CurrentTimeZone.GetUtcOffset([datetime]::Now).TotalHours          (yields integer)
    # Unix: date +%z   =   +hhmm numeric time zone (e.g., -0400)
    UTC_N_hours_offset = as.integer(substr(Sys.getenv("HOST_TIMEZONE_UTC_OFFSET"), 1, 3))

    if(!is.na(UTC_N_hours_offset)) {
      timestamp_prettyprint = format(Sys.time() + UTC_N_hours_offset * 3600, format = "%Y-%m-%d_%H;%M;%S", tz = "UTC")
    } else {
      timestamp_prettyprint = format(Sys.time(), format = "%Y-%m-%d_%H;%M;%S")
    }
    output_dir = path_clean_slashes(paste0(output_dir, "/", timestamp_prettyprint, "/"))
  }

  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
  if(!dir.exists(output_dir)) {
    append_log(paste("failed to create directory;", output_dir), type = "error")
  }

  # check if all output files we want to write to are accessible, before running all time-consuming code
  fname_samples = path_append_and_check(output_dir, "samples.xlsx")
  remove_file_if_exists(fname_samples)

  fname_stats = path_append_and_check(output_dir, "differential_abundance_analysis.xlsx")
  remove_file_if_exists(fname_stats)

  if(output_qc_report) {
    remove_file_if_exists(path_append_and_check(output_dir, "report.pdf"))
  }
  if(output_abundance_tables) {
    fname_abundances = path_append_and_check(output_dir, "peptide_and_protein_abundances.xlsx")
    remove_file_if_exists(fname_abundances)
  }


  # facilitate multiprocessing, this is only used for our custom implementation of msqrob at the moment
  if(any(c("msqrob", "msqrobsum") %in% dea_algorithm)) {
    # speed up your analysis by using multiple processor cores (default; all cores but one)
    cl <<- initialize_multiprocessing()
    # finally, shut down multithreading clusters. use on.exit to execute regardless of downstream errors
    on.exit({ suppressWarnings(parallel::stopCluster(cl)); rm(cl, envir = .GlobalEnv) })
  }


  ##### all filters; global, local, by_group (for CoV plots). adds additional columns to peptide table with filtered and normalized peptide intensities
  dataset = filter_dataset(dataset,
                           filter_min_detect = filter_min_detect, filter_fraction_detect = filter_fraction_detect, filter_min_quant = filter_min_quant, filter_fraction_quant = filter_fraction_quant,
                           filter_min_peptide_per_prot = filter_min_peptide_per_prot, filter_topn_peptides = filter_topn_peptides,
                           norm_algorithm = norm_algorithm,
                           by_group = T, all_group = (!filter_by_contrast || output_abundance_tables || output_qc_report), by_contrast = filter_by_contrast)


  ##### DE analysis
  # quantitative analysis; eBayes/MSqRob/MS-EmpiRe/MSqRobSum
  dataset = dea(dataset, qval_signif = dea_qvalue_threshold, fc_signif = dea_log2foldchange_threshold, algo_de = dea_algorithm, algo_rollup = dea_protein_rollup, output_dir_for_eset = ifelse(dump_all_data, output_dir, ""))
  # debug; print(dataset$de_proteins %>% filter(signif) %>% left_join(dataset$proteins))

  # qualitative analysis
  dataset = differential_detect(dataset, min_samples_observed = diffdetect_min_samples_observed)

  # export data tables to file
  openxlsx::write.xlsx(dataset$samples %>% select(!!grep("^key_", colnames(dataset$samples), ignore.case = T, value = T, invert = T)), fname_samples)
  export_statistical_results(dataset, output_dir)
  if(output_abundance_tables) {
    export_protein_abundance_matrix(dataset, algo_rollup = dea_protein_rollup, output_dir = output_dir)
  }

  # write all data tables to compressed .tsv files
  if(dump_all_data) {
    save(dataset, file = path_append_and_check(output_dir, "dataset.RData"), compress = T)
    data.table::fwrite(dataset$peptides, path_append_and_check(output_dir, "peptides.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$proteins, path_append_and_check(output_dir, "proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$samples, path_append_and_check(output_dir, "samples.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    if(is_tibble(dataset$de_proteins) && nrow(dataset$de_proteins) > 0) {
      data.table::fwrite(dataset$de_proteins, path_append_and_check(output_dir, "de_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    }
    if(is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0) {
      data.table::fwrite(dataset$dd_proteins, path_append_and_check(output_dir, "dd_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    }
  }

  # QC report
  if(output_qc_report) {
    generate_pdf_report(dataset, output_dir = output_dir, pca_sample_labels = pca_sample_labels)
  }

  append_log(paste("output directory;", output_dir), type = "info")
  return(dataset)
}
