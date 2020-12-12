
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
generate_pdf_report = function(dataset, output_dir, norm_algorithm = "vwmb", pca_label_samples_by_shortname = NA) {

  start_time = Sys.time()
  append_log("creating PDF report...", type = "progress")

  ################ prepare data ################
  isdia = is_dia_dataset(dataset)

  # basic filtering within each group
  dataset = filter_peptides_by_group(dataset, colname="intensity_qc_basic", disregard_exclude_samples = FALSE, nquant=2, ndetect=1+isdia, norm_algorithm = norm_algorithm)

  # for DIA, we remove non-detect observations from the tibble used for RT QC plots (poor Q-value may result in large RT/intensity deviation, especially since we here filter very loosely)
  if(isdia) {
    dataset$peptides$intensity_qc_basic[dataset$peptides$detect != TRUE] <- NA
  }

  # we need by-group filtering downstream for CoV plots
  # note; cannot recycle intensity_qc_basic data since this should discard outlier samples and keep MBR intensities (eg; peptide*sample intensity values but not 'detected')
  if(!"intensity_by_group" %in% colnames(dataset$peptides)) {
    append_log(sprintf("generating 'by_group' filtered&normalized data for downstream CoV plots, filter min-detect:%d min-quant:2 ...", 1+isdia), type = "progress")
    dataset = filter_peptides_by_group(dataset, colname="intensity_by_group", disregard_exclude_samples = TRUE, nquant=2, ndetect=1+isdia, norm_algorithm = norm_algorithm)
  }

  # need all-group filtering downstream for PCA plots
  if(!"intensity_all_group" %in% colnames(dataset$peptides)) {
    append_log(sprintf("generating 'all_group' filtered&normalized data for downstream PCA plots, filter min-detect:%d min-quant:2 ...", 0+isdia), type = "progress")
    # inefficient code, but simple / easy to read for now. Rare use-case anyway (because user didn't use standard pipeline nor filter_dataset upstream)
    dataset2 = filter_dataset(dataset, filter_min_detect = 0+isdia, filter_min_quant = 2, norm_algorithm = norm_algorithm, by_group = F, by_contrast = F, all_group = T)
    dataset$peptides$intensity_all_group = dataset2$peptides$intensity_all_group
    rm(dataset2)
  }


  ## optionally, we can use alternative filter for peptides to be used in PCA
  # either recycle user settings, or follow below filter and add fraction_detect and fraction_quant
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
  samples_colors_long = dataset$samples %>% select(sample_id, shortname) %>%
    left_join(sample_color_coding(dataset$samples), by="sample_id") %>%
    # set sample property factor levels (arrangement for plotting later) according to the column name order in the sample metadata table
    mutate(prop = factor(prop, levels=intersect(colnames(dataset$samples), prop))) %>%
    # sort such that first, properties are arranged accordingly to factor levels and then samples by their order in the sample metadata table
    arrange(as.numeric(prop), match(sample_id, dataset$samples$sample_id))

  # for convenience in some plotting functions, the color-codings as a wide-format tibble
  samples_colors = samples_colors_long %>% select(sample_id, shortname, prop, clr) %>% pivot_wider(id_cols = c(sample_id, shortname), names_from=prop, values_from=clr)


  ggplot_cscore_histograms = list()
  if(length(dataset$plots) > 0 && length(dataset$plots$ggplot_cscore_histograms) > 0) {
    ggplot_cscore_histograms = dataset$plots$ggplot_cscore_histograms
  }

  ### contrasts
  append_log("report: constructing plots specific for each contrast", type = "progress")
  l_contrast = list()
  if(is_tibble(dataset$de_proteins) && nrow(dataset$de_proteins) > 0) {
    de_proteins = dataset$de_proteins %>% left_join(dataset$proteins, by="protein_id")

    column_contrasts = dataset_contrasts(dataset)
    for(contr in column_contrasts) {
      stats_contr = de_proteins %>% filter(contrast == contr)

      mtitle = contr
      if("contrast_ranvars" %in% colnames(stats_contr) && length(stats_contr$contrast_ranvars[1]) > 0 && !is.na(stats_contr$contrast_ranvars[1]) && stats_contr$contrast_ranvars[1] != "") {
        mtitle = paste0(mtitle, "\nuser-specified random variables added to regression model: ", stats_contr$contrast_ranvars[1])
      }

      # optionally, provide thresholds for foldchange and qvalue so the volcano plot draws respective lines
      l_contrast[[contr]] = list(p_volcano_contrast = plot_volcano(stats_de = stats_contr, log2foldchange_threshold = ifelse(stats_contr$signif_threshold_log2fc[1]==0, NA, stats_contr$signif_threshold_log2fc[1]), qvalue_threshold = stats_contr$signif_threshold_qvalue[1], mtitle = mtitle),
                                 p_pvalue_hist = plot_pvalue_histogram(stats_contr %>% mutate(color_code = algo_de), mtitle=contr),
                                 p_foldchange_density = plot_foldchanges(stats_contr %>% mutate(color_code = algo_de), mtitle=contr))
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
      arrange(match(contrast, column_contrasts), algorithm) %>%
      ungroup() %>%
      mutate(contrast = sub("^contrast: ", "", contrast))

    # optionally, limit contrast string length (evenly on each side by N characters)
    tib_report_stats_summary$contrast = unlist(lapply(strsplit(tib_report_stats_summary$contrast, " vs ", fixed = T),
                                                      function(x) paste(stringr::str_trunc(x, 18, "right"), collapse = " vs ")))


    ### analogous for differential detection
    if(is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0) {
      tib_report_diffdetects_summary = dataset$dd_proteins %>%
        filter(!is.na(zscore_count_detect)) %>%
        left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by="protein_id") %>%
        # sort such that top hits come first
        arrange(desc(abs(zscore_count_detect))) %>%
        # summary stats per contrast
        group_by(contrast) %>%
        summarise(`#proteins tested` = n(),
                  `#abs(zscore) >= 3` = sum(abs(zscore_count_detect) >= 3),
                  `top10` = tolower(paste(stringr::str_trunc(head(gene_symbols_or_id, 10), width = 10, side = "right"), collapse=", ") )) %>%
        ungroup() %>%
        arrange(match(contrast, column_contrasts)) %>%
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
      # downstream we use a helper function to (naively) remove comments and compactly format the code
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

  # rmarkdown known bugs: do not use output_dir  @  https://github.com/rstudio/rmarkdown/issues/861
  # a recommended work-around is to set the base/root directories inside the rmarkdown report using knitr::opts_chunk$set(...)
  # ref; https://github.com/yihui/knitr/issues/277#issuecomment-6528846
  # our robust work-around: since the render function uses the report document and it's dir as a base, we simply copy the Rmarkdown file to output dir and run render() with all default settings
  output_dir__temp = paste0(output_dir, "/temp", floor(as.numeric(Sys.time())) ) # to ensure unique dirname, simply add unix timestamp with seconds as precision
  dir.create(output_dir__temp, showWarnings = F) # don't care if directory already exists
  if(!dir.exists(output_dir__temp)) {
    append_log(paste("failed to create temp directory at;", output_dir__temp), type = "error")
  }

  f_newlocation = paste0(output_dir__temp, "/", basename(f))
  fpdf_newlocation = paste0(output_dir__temp, "/report.pdf")
  # copy .Rmd file into temp directory, nested in chosen output dir
  if(!file.copy(from = f, to = f_newlocation)) {
    append_log(paste("failed to move report into to the output directory:", fpdf_newlocation), type = "error")
  }

  # create the actual report
  rmarkdown::render(input = f_newlocation, output_file = "report.pdf", quiet = T, clean = TRUE) # , output_file = paste0(output_dir, "/report.pdf") , intermediates_dir = output_dir, knit_root_dir = output_dir

  # move report to output dir and remove temp dir
  if(!file.exists(fpdf_newlocation)) {
    append_log(paste("there is no report at:", fpdf_newlocation), type = "error")
  }

  file.rename(fpdf_newlocation, paste0(output_dir, "/report.pdf"))
  # try to remove our temp files
  file.remove(f_newlocation)
  file.remove(paste0(output_dir__temp, "/report.tex"))
  # try to remove entire temp dir; may fail if user opened one of the files or is inside the dir in windows explorer
  # should be safe because we use a unique name in a dir we created previously AND we checked that this is an existing path where we have write access (otherwise above code would have failed)
  unlink(output_dir__temp, recursive = T, force = T) # use recursive=T, because unlink() documentation states: "If recursive = FALSE directories are not deleted, not even empty ones"

  append_log_timestamp("report:", start_time)
}
