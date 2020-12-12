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
#' all-in-one function that covers the vast majority of use-cases. Typically, this function is all you need (besides importing your data upstream).
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
#' filter_min_detect = 1, filter_fraction_detect = 0.25, filter_min_quant = 3, filter_fraction_quant = 0.75
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
#' options: "" (empty string disables normalization), "vsn", "loess", "rlr", "msempire", "vwmb", "modebetween".
#'
#' You can combine normalizations by providing an array of options to apply subsequential normalizations.
#'
#' For instance, \code{norm_algorithm = c("vsn", "modebetween")} applies the vsn algorithm (quite strong normalization reducing variation) and then balances between-group foldchanges with the modebetween normalization from the vwmb algorithm.
#'
#' By default uses the built-in normalization algorithm Variation Within, Mode Between (vwmb).
#' Benchmarks have shown that vsn&modebetween (code example above) is a great alternative to consider.
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
#' options: ebayes, msempire, msqrob (to run multiple, provide an array)
#'
#' Users are encouraged to simply run all statistical models and consider the set of candidate proteins that pop up in multiple approaches. Output tables already contain columns to convenience this approach, if 3 DEA algorithms have been selected.
#' e.g.; \code{dea_algorithm = c("ebayes", "msempire", "msqrob")}
#'
#' All statistical results are shown in the QC report and output tables and can be easily compared.
#'
#'
#'
#' @param dataset a valid dataset object generated upstream by, for instance, import_dataset_skyline
#' @param filter_min_detect in order for a peptide to 'pass' in a sample group, in how many replicates must it be detected?
#' @param filter_fraction_detect in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be detected?
#' @param filter_min_quant in order for a peptide to 'pass' in a sample group, in how many replicates must it be quantified?
#' @param filter_fraction_quant in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be quantified?
#' @param filter_min_peptide_per_prot in order for a peptide to 'pass' in a sample group, how many peptides should be available after detect filters?
#' @param filter_topn_peptides maximum number of peptides to maintain for each protein (from the subset that passes above filters, peptides are ranked by the number of samples where detected and their variation between replicates). 1 is default, 2 can be a good choice situationally. If set to 1, make sure to inspect individual peptide data/plots for proteins with 1 peptide.
#' @param filter_by_contrast should the above filters be applied to all sample groups, or only those tested within each contrast? Enabling this optimizes available data in each contrast, but increases the complexity somewhat as different subsets of peptides are used in each contrast and normalization is applied separately.
#' @param norm_algorithm normalization algorithms. Can be empty string to skip normalization, or any algorithm listed by msdap::normalization_algorithms(). Provide an array of options to run each algorithm consecutively, for instance; c("vsn", "modebetween_protein") to first apply vsn normalization and then correct between-group ratios such that the protein-level log2-foldchange mode is zero
#' @param dea_protein_rollup strategy for combining peptides to proteins as used in DEA algorithms that first combine peptides to proteins and then apply statistics, like ebayes and deqms. options: maxlfq, sum. The former applies the MaxLFQ algorithm, the latter employs the 'classic' strategy of summing all peptides per protein
#' @param dea_algorithm algorithms for differential expression analysis. options: ebayes, deqms, msqrobsum, msempire, msqrob (to run multiple, provide an array)
#' @param dea_qvalue_threshold threshold for significance of adjusted p-values in figures and output tables
#' @param dea_log2foldchange_threshold threshold for significance of log2 foldchanges. Set to zero to disregard, a positive value to apply a cutoff to absolute foldchanges or use bootstrap analyses to infer a suitable foldchange threshold by providing either NA or a negative value. default: 0
#' @param diffdetect_min_samples_observed for differential detection only; minimum number of samples where a protein should be observed at least once by any of its peptides (in either group) when comparing a contrast of group A vs B
#' @param plot_pca_label_by_shortname whether to use sample names or a numeric ID as labels in the PCA plot. options: NA (let the code decide, default), TRUE (always use sample 'shortname'), FALSE (always use numeric ID)
#' @param output_dir output directory where all output files are stored, must be an existing directory
#' @param output_within_timestamped_subdirectory optionally, automatically create a subdirectory (without output_dir) that has the current date&time as name and store results there
#' @param output_abundance_tables whether to create an Excel document with all peptide abundances, with multiple sheets indicating results from all filters applied. For large datasets this results in huge files, so only recommended for small datasets. options: FALSE, TRUE
#' @param output_qc_report whether to create the Quality Control report. options: FALSE, TRUE . Highly recommended to set to TRUE. Set to FALSE to skip the report PDF (eg; to just do statistics and skip the time-consuming report creation)
# @param output_peptide_plots whether to create a plot for each protein detailing all of it's peptide abundances. options: "" or "none" to disable, "signif"=significant proteins only, "all"=all proteins, "complete"=all proteins and split by DE algorithm
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
                               # normalization
                               norm_algorithm = c("vwmb", "modebetween_protein"),
                               # DE algorithms
                               dea_protein_rollup = "maxlfq",
                               dea_algorithm = c("deqms", "msqrob"), # options; ebayes, deqms, msqrobsum, msempire, msqrob (to run multiple, provide an array)
                               # significance cutoff
                               dea_qvalue_threshold = 0.01,
                               dea_log2foldchange_threshold = 0,
                               # differential detection
                               diffdetect_min_samples_observed = 3,
                               # plot options
                               plot_pca_label_by_shortname = NA, # if NA, auto estimate
                               #
                               output_abundance_tables = FALSE,
                               output_qc_report = TRUE,
                               # output_peptide_plots = "none", # options: "signif"=significant proteins only, "all"=all proteins, "complete"=all proteins and split by DE algorithm
                               # where to store results
                               output_dir,
                               output_within_timestamped_subdirectory = TRUE,
                               dump_all_data = FALSE) {

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
    output_dir = paste0(output_dir, "/", timestamp_prettyprint, "/")
  }

  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
  if(!dir.exists(output_dir)) {
    append_log(paste("failed to create directory;", output_dir), type = "error")
  }

  # check if all output files we want to write to are accessible, before running all time-consuming code
  fname_samples = sprintf("%s/samples.xlsx", output_dir)
  fname_abundances = sprintf("%s/peptide_and_protein_abundances.xlsx", output_dir)

  remove_file_if_exists(fname_samples)
  remove_file_if_exists(paste0(output_dir, "/differential_abundance_analysis.xlsx"))

  if(output_qc_report) {
    remove_file_if_exists(paste0(output_dir, "/report.pdf"))
  }
  if(output_abundance_tables) {
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

  # write analysis data tables to file
  openxlsx::write.xlsx(dataset$samples %>% select(!!grep("^key_", colnames(dataset$samples), ignore.case = T, value = T, invert = T)), fname_samples)
  write_statistical_results_to_file(dataset, output_dir)
  if(output_abundance_tables) {
    # only do mode normalization, no normalizations that 'transform' (eg; vsn), users should do do downstream if they want to filter etc first
    write_peptide_abundance_matrix_to_file(dataset, filename = fname_abundances, norm_algorithm = norm_algorithm)
  }

  # write all data tables to compressed .tsv files
  if(dump_all_data) {
    save(dataset, file = paste0(output_dir, "/dataset.RData"), compress = T)
    data.table::fwrite(dataset$peptides, paste0(output_dir, "/peptides.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$proteins, paste0(output_dir, "/proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$samples, paste0(output_dir, "/samples.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$de_proteins, paste0(output_dir, "/de_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$dd_proteins, paste0(output_dir, "/dd_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
  }

  # QC report
  if(output_qc_report) {
    generate_pdf_report(dataset, output_dir = output_dir, pca_label_samples_by_shortname = plot_pca_label_by_shortname)
  }

  # # peptide data plots; "signif" | "all" | "complete"
  # if(length(output_peptide_plots) == 0 || (length(output_peptide_plots) == 1 && (is.na(output_peptide_plots) || output_peptide_plots == "" || output_peptide_plots == "none"))) {
  #   # valid option, don't plot
  # } else {
  #   # we expect valid user parameter at this point
  #   if(length(output_peptide_plots) == 1 && is.character(output_peptide_plots) && output_peptide_plots %in% c("signif", "all", "complete")) {
  #     plot_peptide_data(dataset$peptides, dataset$proteins, dataset$samples, dataset$de_proteins, output_dir, plot_each_statistical_approach = (output_peptide_plots == "complete"), plot_all_proteins = (output_peptide_plots %in% c("all", "complete")), plot_comparison_msqrob_msempire = FALSE)
  #   } else {
  #     append_log(paste("invalid value(s) for parameter 'output_peptide_plots' (expected any of; NA, 'none', 'signif', 'all', 'complete') but got:", paste(output_peptide_plots, collapse=",")), type = "warning")
  #   }
  # }

  return(dataset)
}
