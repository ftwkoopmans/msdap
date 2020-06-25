
#' Create the Quality Control report
#'
#' Normally used as part of the analysis_quickstart() function, refer to its implementation for example code.
#'
#' @param dataset a valid dataset object generated upstream by, for instance, import_dataset_skyline
#' @param output_dir output directory where all output files are stored, must be an existing directory
#' @param norm_algorithm normalization algorithm(s) used to normalize the QC proportion of data. options; "", "vsn", "loess", "rlr", "msempire", "vwmb", "modebetween". Provide an array of options to run each algorithm consecutively
#' @param pca_label_samples_by_shortname whether to use sample names or a numeric ID as labels in the PCA plot. options: NA (let the code decide, default), TRUE (always use sample 'shortname'), FALSE (always use numeric ID)
#'
#' @import knitr
#' @importFrom rmarkdown render
#' @importFrom devtools session_info
#' @importFrom xtable xtable
#' @importFrom stringr str_wrap
#' @importFrom devtools session_info
#' @export
generate_pdf_report = function(dataset, output_dir, norm_algorithm = "vwmb", pca_label_samples_by_shortname = TRUE) {

  start_time <- Sys.time()
  append_log("creating PDF report...", type = "progress")

  ################ prepare data ################
  isdia = is_dia_dataset(dataset)

  # basic filtering within each group
  dataset = filter_peptides_by_group(dataset, colname="intensity_qc_basic", disregard_exclude_samples=F, nquant=2, ndetect=1+isdia, norm_algorithm = norm_algorithm)

  # TODO: we still want this? or decide within respective plot functions ?
  # for DIA, we remove non-detect observations from the tibble used for RT QC plots (poor Q-value may result in large RT/intensity deviation, especially since we here filter very loosely)
  if(isdia) {
    dataset$peptides$intensity_qc_basic[dataset$peptides$detect != TRUE] <- NA
  }


  # # for CoV leave-one-out strategy, filter by-group following user settings. Note that outlier samples are included, just like with within-group peptide intensity fold-change analysis !
  # if(any(table(dataset$samples$group) >= 4)) { # skip if dataset has no groups with minimum number of required replicates
  #   dataset = filter_peptides_by_group(dataset, colname="intensity_qc_cov_loo", disregard_exclude_samples=F, nquant=3, ndetect=ifelse(isdia, 3, 1), norm_algorithm = norm_algorithm)
  # }

  # TODO: settle on final PCA filter. either recycle user settings, or follow below filter and add fraction_detect and fraction_quant
  # !! these columns are not used downstream @ report2.Rmd yet
  # if(any(dataset$samples$exclude) && "intensity_all_group_withexclude" %in% colnames(dataset$peptides) && any(!is.na(dataset$peptides$intensity_all_group_withexclude))) {
  #   # there are exclude samples, do something with this column...
  # }
  # # dataset = filter_peptides_by_group(dataset, colname="intensity_qc_pca_all", disregard_exclude_samples=F, nquant=3, ndetect=ifelse(isdia, 3, 1), norm_algorithm = "vwmb") # intensity_all_group_withexclude
  # # dataset = filter_peptides_by_group(dataset, colname="intensity_qc_pca_noexclude", disregard_exclude_samples=T, nquant=3, ndetect=ifelse(isdia, 3, 1), norm_algorithm = "vwmb") # intensity_all_group
  # dataset$peptides$intensity_qc_pca_all = dataset$peptides$intensity_all_group_withexclude
  # dataset$peptides$intensity_qc_pca_noexclude = dataset$peptides$intensity_all_group



  ################ plot ################

  ### sample metadata, color-coding and generic ggplots
  samples_colors = sample_color_coding(dataset$samples)
  samples_colors_long = sample_color_coding_as_long_tibble(dataset$samples, samples_colors)

  ggplot_cscore_histograms = list()
  if(length(dataset$plots) > 0 && length(dataset$plots$ggplot_cscore_histograms) > 0) {
    ggplot_cscore_histograms = dataset$plots$ggplot_cscore_histograms
  }

  ### contrasts
  append_log("report: constructing plots specific for each contrast", type = "progress")
  l_contrast = list()
  if(is_tibble(dataset$de_proteins) && nrow(dataset$de_proteins) > 0) {
    de_proteins = dataset$de_proteins %>% left_join(dataset$proteins, by="protein_id")

    column_contrasts = grep("^contrast:", colnames(dataset$samples), ignore.case = T, value = T)
    for(contr in column_contrasts) {
      stats_contr = de_proteins %>% filter(contrast == contr)
      # optionally, provide thresholds for foldchange and qvalue so the volcano plot draws respective lines
      l_contrast[[contr]] = list(p_volcano_contrast = suppressWarnings(plot_volcano(stats_de = stats_contr, log2foldchange_threshold = ifelse(stats_contr$signif_threshold_log2fc[1]==0, NA, stats_contr$signif_threshold_log2fc[1]), qvalue_threshold = stats_contr$signif_threshold_qvalue[1], mtitle = contr)),
                                 p_pvalue_hist = suppressWarnings(plot_pvalue_histogram(stats_contr %>% mutate(color_code = algo_de), mtitle=contr)),
                                 p_foldchange_density = suppressWarnings(plot_foldchanges(stats_contr %>% mutate(color_code = algo_de), mtitle=contr)))
    }


    ### summary table of all stats
    tib_report_stats_summary = de_proteins %>%
      arrange(pvalue) %>%
      rename(algorithm = algo_de) %>%
      # filter(signif) %>% # include this to remove contrast*algo_de entries without any significant hits
      group_by(contrast, algorithm) %>%
      summarise(`#test` = sum(is.finite(pvalue)),
                `#hits` = sum(signif),
                `top10 significant` = tolower(paste(head(stringr::str_trunc(gene_symbols_or_id[signif], 10, "right"), 10), collapse=", ") )) %>%
      # summarise(`#tested` = sum(is.finite(pvalue)), `#signif` = sum(signif),
      #           `FC>=1.1` = sum(signif & abs(foldchange.log2) >= log2(1.1), na.rm=T), `FC>=1.5` = sum(signif & abs(foldchange.log2) >= log2(1.5), na.rm=T) ) %>%
      arrange(contrast, algorithm) %>%
      ungroup() %>%
      mutate(contrast = sub("^contrast: ", "", contrast))

    # optionally, limit contrast string length (evenly on each side by N characters)
    tib_report_stats_summary$contrast = unlist(lapply(strsplit(tib_report_stats_summary$contrast, " vs ", fixed = T),
                                                      function(x) paste(stringr::str_trunc(x, 18, "right"), collapse = " vs ")))


    ### analogous for differential detection
    if(is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0) {
      tib_report_diffdetects_summary = dataset$dd_proteins %>%
        left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by="protein_id") %>%
        filter(!is.na(diff_detect_zscore)) %>%
        arrange(desc(abs(diff_detect_zscore))) %>%
        group_by(contrast) %>%
        summarise(`#test` = n(),
                  `#hits` = sum(diff_detect_zscore_candidate),
                  `top10 candidates` = tolower(paste(head(stringr::str_trunc(gene_symbols_or_id[diff_detect_zscore_candidate], 10, "right"), 10), collapse=", ") )) %>%
        arrange(contrast) %>%
        ungroup() %>%
        mutate(contrast = sub("^contrast: ", "", contrast))

      # optionally, limit contrast string length (evenly on each side by N characters)
      tib_report_diffdetects_summary$contrast = unlist(lapply(strsplit(tib_report_diffdetects_summary$contrast, " vs ", fixed = T),
                                                              function(x) paste(stringr::str_trunc(x, 18, "right"), collapse = " vs ")))
    }
  }


  ### differential detect
  if("dd_proteins" %in% names(dataset) && is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0) {
    dd_plots = plot_differential_detect(dataset)
  }



  ################ history ################
  history_as_string = ""
  if(exists("history_")) {
    history_as_string = history_
  } else {
    try({
      ## save command history as variable, so we can show in report
      f_hist = paste0(output_dir, "/.Rhistory")
      # !! obscure bug; if we use utils::savehistory() the function silently fails, usage without the package prefix work fine. importing utils in namespace also seems to give problems
      savehistory(f_hist)
      # use a helper function to (naively) remove comments and compactly format the code
      history_as_string = tryCatch(readLines(f_hist), warning = function(x){""}, error = function(x){""})
      if(file.exists(f_hist)) {
        file.remove(f_hist)
      }
    }, silent = TRUE)
  }


  ################ render Rmarkdown ################

  f = system.file("rmd", "report.Rmd", package = "msdap")
  if(!file.exists(f)) {
    append_log(paste("cannot find report template file:", f), type = "error")
  }
  append_log("report: rendering report (this may take a while depending on dataset size)", type = "progress")
  rmarkdown::render(f, output_dir = output_dir, output_file = "report.pdf", quiet = T, intermediates_dir = output_dir, knit_root_dir = output_dir)
  append_log_timestamp("report:", start_time)
}
