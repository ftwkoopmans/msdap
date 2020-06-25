
#' iterate all contrasts: gather peptide data used in upstream stats, create individual protein plots, finally merge plots into end result PDFs
#' temp filenames not random, so only run a single session of this tool at the same time on a given system
#'
#' @param peptides todo
#' @param proteins proteins tibble needs columns: (protein_id, fasta_headers)
#' @param samples todo
#' @param de_proteins todo
#' @param output_dir todo
#' @param plot_each_statistical_approach todo
#' @param plot_all_proteins todo
#' @param plot_comparison_msqrob_msempire todo
plot_peptide_data = function(peptides, proteins, samples, de_proteins, output_dir, plot_each_statistical_approach = TRUE, plot_all_proteins = FALSE, plot_comparison_msqrob_msempire = FALSE) {

  # peptides; sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, detect, intensity=!!col_contr_intensity
  # proteins; protein_id, fasta_headers
  # samples; sample_id, shortname, group, condition = !!col_contr
  # de_proteins; "pvalue" "qvalue" "signif" "foldchange.log2" "protein_id" "algo_de" "contrast"

  if(!is_tibble(de_proteins) || nrow(de_proteins) == 0) {
    append_log("no differential expression analysis results, not plotting peptide-level data", type = "warning")
    return()
  }

  column_contrasts = grep("^contrast:", colnames(samples), ignore.case = T, value = T)
  append_log("plot peptide data for each protein...", type = "info")
  for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
    start_time <- Sys.time()
    lbl = sub("^contrast: *", "", col_contr)

    samples_for_contrast = samples
    ## enable this filter step to plot only samples that match current contrast
    # samples_for_contrast = samples %>%
    #   select(sample_id, shortname, group, condition = !!col_contr) %>%
    #   filter(condition != 0) %>%
    #   arrange(condition)

    col_contr_intensity = get_column_intensity(peptides, col_contr)
    append_log(paste("using data from peptide filter:", names(col_contr_intensity)), type = "info")
    tib_contrast_peptides = peptides %>%
      select(sample_id, protein_id, peptide_id, sequence_plain, sequence_modified, detect, intensity=!!as.character(col_contr_intensity)) %>%
      filter(sample_id %in% samples_for_contrast$sample_id & is.finite(intensity))

    tib_contrast_stats = de_proteins %>% filter(contrast == col_contr) %>% arrange(pvalue)

    # select proteins to visualize; if user wants all proteins, select everything, otherwise only significant
    protein_id_plot = tib_contrast_stats %>% filter(plot_all_proteins | signif) %>% distinct(protein_id) %>% pull()

    if(length(protein_id_plot) == 0)
      next

    # plot each protein, immediately storing each ggplot in a pdf file (RAM efficient + enables recycling of plots by simply merging individual PDFs downstream)
    # result tibble contains the protein_id's and their respective pdf filenames
    tib_protein_pdf_files = plot_peptide_data_per_protein(peptides = tib_contrast_peptides, proteins = proteins, stats_de = tib_contrast_stats,
                                                          samples = samples, protein_ids = protein_id_plot)

    ## assemble various PDFs with specific plots-of-interest. eg; all proteins significant by any method, detected by method X  or  by method X and not method Y
    if(plot_all_proteins) {
      pdf_combine_chunks(input = tib_protein_pdf_files$pdf_filename, output = sprintf("%s/all_proteins - %s.pdf", output_dir, lbl))
    }
    pdf_combine_chunks(input = tib_protein_pdf_files %>% filter(protein_id %in% (tib_contrast_stats %>% filter(signif) %>% pull(protein_id)) ) %>% pull(pdf_filename),
                output = sprintf("%s/all_significant_proteins - %s.pdf", output_dir, lbl))

    # list of significant proteins per DE algorithm (already sorted by pvalue)
    ualg = unique(tib_contrast_stats$algo_de)
    if(plot_each_statistical_approach && length(ualg) > 1) {
      for(alg in ualg) {
        alg_protein_id = tib_contrast_stats %>% filter(algo_de == alg & signif) %>% pull(protein_id)
        if(length(alg_protein_id) > 0) {
          pdf_combine_chunks(input = tib_protein_pdf_files %>% filter(protein_id %in% alg_protein_id) %>% pull(pdf_filename),
                             output = sprintf("%s/%s_significant - %s.pdf", output_dir, alg, lbl))
        }
      }
    }

    # unique hits in msempire vs msqrob, and vice versa
    if(plot_comparison_msqrob_msempire && all(c("msempire", "msqrob") %in% ualg)) {
      msempire_protein_id = tib_contrast_stats %>% filter(algo_de == "msempire" & signif) %>% pull(protein_id)
      msqrob_protein_id = tib_contrast_stats %>% filter(algo_de == "msqrob" & signif) %>% pull(protein_id)
      msempire_pdf = tib_protein_pdf_files %>% filter(protein_id %in% setdiff(msempire_protein_id, msqrob_protein_id)) %>% pull(pdf_filename)
      msqrob_pdf = tib_protein_pdf_files %>% filter(protein_id %in% setdiff(msqrob_protein_id, msempire_protein_id)) %>% pull(pdf_filename)
      if(length(msempire_pdf) > 0) {
        pdf_combine_chunks(input = msempire_pdf, output = sprintf("%s/msempire_significant_not_msqrob - %s.pdf", output_dir, lbl))
      }
      if(length(msqrob_pdf) > 0) {
        pdf_combine_chunks(input = msqrob_pdf, output = sprintf("%s/msqrob_significant_not_msempire - %s.pdf", output_dir, lbl))
      }
    }

    append_log_timestamp(paste(col_contr, ", plot peptide data per protein"), start_time)
  }
}


# TODO
# # implementation analogous to plot_proteins_significant()
# plot_proteins_of_interest = function(stats, samples, proteins, output_file, protein_ids) {
#   column_contrasts = grep("^contrast:", colnames(samples), ignore.case = T, value = T)
#   for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
#     lbl = sub("^contrast: *", "", col_contr)
#     tib_contrast_peptides = stats$peptides %>% filter(contrast == col_contr) %>% select(-intensity) %>% rename(intensity=intensity_norm)
#     tib_contrast_stats = stats$stats %>% filter(contrast == col_contr) %>% arrange(pvalue)
#
#     if(any(protein_ids %in% tib_contrast_stats$protein_id)) {
#       tib_protein_pdf_files = plot_peptide_data_per_protein(peptides = tib_contrast_peptides, proteins = proteins, stats_de = tib_contrast_stats,
#                                                             samples = samples, protein_ids = protein_ids)
#       pdf_combine_chunks(input = tib_protein_pdf_files$pdf_filename, output = sprintf("%s - %s.pdf", remove_file_extension_from_path(output_file), lbl))
#     }
#   }
# }
#



#' placeholder title
#' proteins tibble needs columns: (protein_id, fasta_headers)
#' @param peptides todo
#' @param proteins todo
#' @param stats_de todo
#' @param samples todo
#' @param protein_ids todo
#'
#' @importFrom foreach foreach %dopar%
plot_peptide_data_per_protein = function(peptides, proteins, stats_de, samples, protein_ids) {
  # append_log("plotting peptide data and stats per protein...", type = "info")

  ### prep input data

  # subset user-requested proteins to those present in data
  protein_ids = intersect(protein_ids, peptides$protein_id)
  tib_proteins = proteins %>%
    select(protein_id, fasta_headers) %>%
    filter(protein_id %in% protein_ids) %>%
    mutate(pid=protein_id) %>%
    arrange(match(protein_id, protein_ids)) %>%
    # add_column(filename = tempfile(paste0("ppipeline_prot_", seq_along(protein_ids))))
    add_column(filename = sprintf("%s/ppipeline_prot_%s.pdf", tempdir(), seq_along(protein_ids))) # non-unique, but we don't care as we encourage over-writing by the same tool to reduce disk space. TODO: should we include some unique session-id or cleanup files downstream

  if(!all(tib_proteins$protein_id == protein_ids)) {
    append_log("all protein ids have to be present in protein metadata table", type = "error")
  }

  tib_sample_groups = samples %>% select(sample_id, shortname, group) %>% filter(sample_id %in% peptides$sample_id)
  # print(tib_sample_groups)

  # use select() to minimize data columns added to tib_input, which (slightly) reduces the RAM footprint
  tib_input = left_join(peptides %>% select(peptide_id, protein_id, sample_id, detect, intensity) %>% filter(protein_id %in% protein_ids),
                        tib_sample_groups, by = "sample_id")
  # tib_input = left_join(tib_input, proteins %>% select(protein_id, fasta_headers), by = "protein_id")
  tib_input$shortname = factor(tib_input$shortname, levels = tib_sample_groups$shortname)
  tib_input$detect = factor(as.character(tib_input$detect), levels = c("TRUE", "FALSE"))
  # transform peptide intensities from log2 to log10
  tib_input$intensity = tib_input$intensity / log2(10)

  # add a grouping var for each peptide per group
  tib_input = tib_input %>% unite(peptide_id_group, peptide_id, group, remove = FALSE)

  # stats in pretty-print format
  stats_de_prettyprint = stats_de %>%
    filter(protein_id %in% protein_ids) %>%
    select(protein_id, algorithm = algo_de, foldchange = foldchange.log2, pvalue, qvalue, signif) %>%
    mutate(foldchange = sprintf("%.2f", 2^foldchange),
           pvalue = sprintf("%.2g", pvalue),
           qvalue = sprintf("%.2g", qvalue) )


  ### list object with all data required for a plot, per protein
  tib_input = tib_proteins %>% group_by(protein_id) %>%
    inner_join(tib_input %>% group_by(protein_id) %>% nest() %>% rename(peptide_data=data), by = "protein_id") %>%
    left_join(stats_de_prettyprint %>% group_by(protein_id) %>% nest() %>% rename(stat_data=data), by="protein_id") %>%
    nest()
  # tib_input = tib_input %>% group_by(protein_id) %>% nest() %>% rename(peptide_data=data) %>%
  #   left_join(stats_de_prettyprint %>% group_by(protein_id) %>% nest() %>% rename(stat_data=data), by="protein_id") %>%
  #   left_join(tib_proteins, by = "protein_id") %>%
  #   nest()

  ### multiprocessing the actual peptide plots
  tib_input$pdf_filename = foreach::foreach(d = tib_input$data, .combine = 'c', .export = "ggplot_peptide", .packages = c("dplyr", "ggplot2")) %dopar% {
    ggplot_peptide(protein_id=d$pid, filename=d$filename, plot_title = d$fasta_headers, tib_peptides = d$peptide_data[[1]], tib_stats = d$stat_data[[1]], tib_sample_groups=tib_sample_groups)
  }

  tib_input %>% select(protein_id, pdf_filename)
}



#' placeholder title
#' @param protein_id todo
#' @param filename todo
#' @param plot_title todo
#' @param tib_peptides todo
#' @param tib_stats todo
#' @param tib_sample_groups todo
#'
#' @importFrom ggpubr ggtexttable ggarrange
ggplot_peptide = function(protein_id, filename, plot_title, tib_peptides, tib_stats, tib_sample_groups) {
  clr_palette = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999") # hardcoded @ brewer.pal(9, "Set1")

  # ! drop factor levels that are not in the data
  tib_peptides$detect = droplevels(tib_peptides$detect)

  ## data for factor-to-geom_point() decoration, for plot legend
  # # named array maps from `detect` to a shape number (1=circle no fill, 16=circle with fill). use only elements that are actually in the data
  # detect_shape = c("TRUE" = 16, "FALSE" = 1)[levels(tib_peptides$detect)]

  upepid = unique(tib_peptides$peptide_id)
  tib_upep = tibble(peptide_id = upepid, label = gsub("(^_)|(#.*)|(\\.)", "", upepid))
  if(length(upepid) > 9) {
    tib_upep$color = colorRampPalette(clr_palette)(length(upepid))
  } else {
    tib_upep$color = clr_palette[seq_along(upepid)]
  }

  ## add mean abundance. left-join ensures we have an entry for each sample per group
  tib_summ_int = tib_peptides %>%
    group_by(group, peptide_id, peptide_id_group) %>%
    summarise(int_avg = mean(intensity, na.rm = T)) %>%
    left_join(tib_sample_groups, by = "group")


  p = ggplot(tib_peptides, aes(x = shortname, y = intensity, group = peptide_id_group, shape = peptide_id, colour = peptide_id)) +
    geom_point(size = 3, show.legend = F) +
    geom_point(size = 3, aes(alpha = detect, fill = peptide_id)) + # note that we hardcode the levels @ detect, so we know in which order the elements are @ guide_legend()
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = array(tib_upep$color, dimnames=list(tib_upep$peptide_id)), labels = tib_upep$label, name="peptide", aesthetics = c("colour", "fill")) +
    scale_shape_manual(values = array(rep(c(22:25), length.out = length(upepid)), dimnames=list(tib_upep$peptide_id)), labels = tib_upep$label, drop = FALSE, name="peptide") +
    scale_alpha_manual(values = c("TRUE" = .8, "FALSE" = 0), name = "detected?", drop = FALSE) +
    # alpha=F  disables alpha legend to save space
    guides(alpha = F, fill = guide_legend(ncol = 2), colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
    labs(x = "", y = "log10 intensity", fill = "", colour = "", shape = "", title = plot_title) + # set same label/name for all properties of 'sequence_id' to rename the respective legend title
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0, size = 7),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = ifelse(length(upepid) < 10, 7, 5))
    )

  p = p + geom_line(data = tib_summ_int, aes(x = shortname, y = int_avg, group = peptide_id_group, colour = peptide_id), show.legend = F) # we need to set `group` in ggplot to connect the line within peptide * sample group

  # if more than 15 peptides, don't show peptide sequences in legend (gonna be a mess)
  if(length(upepid) > 15) {
    p = p + theme(legend.position = "none")
  }

  # add stats
  p_texttable = ggpubr::ggtexttable(as.data.frame(tib_stats), rows = NULL, theme = ggpubr::ttheme("classic"))

  p_out = ggpubr::ggarrange(p, p_texttable, ncol = 1, nrow = 2, heights = c(.75, .25))
  ggplot2::ggsave(file=filename, plot=p_out, width=7, height=9)

  return(filename)
}



#' placeholder title
#' the function pdf_combine() from pdftools/qpdf package has issues (on windows) when combining hundreds of plots. error: "Too many open files"
#' so we here create a wrapper where we combine source files in batches
#' some references; https://stackoverflow.com/questions/40810704/error-with-r-function-download-file-too-many-open-files
#' @param input todo
#' @param output todo
#'
#' @importFrom pdftools pdf_combine
pdf_combine_chunks = function(input, output) {
  # chunks of at most 200 files  &  respective temporary filenames (re-using first input file to get unique file / for unique purpose)
  input_chunks = split(input, ceiling(seq_along(input) / 200))
  output_temp = paste0(input[1], ".chunk", seq_along(input_chunks))
  # combine files within each chunk
  for(i in seq_along(input_chunks)) {
    pdftools::pdf_combine(input_chunks[[i]], output_temp[i])
  }
  # finally, combine chunks and remove temp files
  pdftools::pdf_combine(output_temp, output)
  res = file.remove(output_temp)
}

