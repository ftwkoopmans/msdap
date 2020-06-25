
#' placeholder title
#' @param tib_input todo
#'
#' @importFrom ggpubr ggarrange theme_classic2
plot_dia_cscore_histograms = function(tib_input) {
  if(length(tib_input) < 2 || !is_tibble(tib_input)) {
    return(list())
  }
  if(!all(c("cscore", "isdecoy") %in% colnames(tib_input)) || !any(tib_input$isdecoy)) {
    append_log("to plot target/decoy score distributions, include columns 'cscore' and 'isdecoy' (with at least some decoy entries) in the peptide report", type = "info")
    return(list())
  }

  plotlist = list()
  plotscores = c()

  for (sid in unique(tib_input$sample_id)) {
    tib = tib_input %>% filter(sample_id == sid & is.finite(cscore))
    if(nrow(tib) == 0) {
      append_log(paste("zero peptides with a cscore in this sample;", sid), type = "warning")
      next
    }

    # classify target and decoy
    tib$type = ifelse(is.finite(tib$isdecoy) & tib$isdecoy == FALSE, "target", "decoy")
    n = table(tib$type)
    if(!"target" %in% names(n)) {
      append_log(paste("no target peptides in sample;", sid), type = "warning")
    }
    if(!"decoy" %in% names(n)) {
      append_log(paste("no decoys peptides in sample;", sid), type = "warning")
    }
    # plot limits
    tib_xlim = c(
      min(0.2, quantile(tib$cscore[tib$type == "decoy"], probs = .01, na.rm = T)),
      max(0.8, quantile(tib$cscore[tib$type == "target"], probs = .99, na.rm = T))
    )
    if(is.na(tib_xlim[1])) {
      tib_xlim[1] = min(0.2, quantile(tib$cscore, probs = .01, na.rm = T))
    }
    if(is.na(tib_xlim[2])) {
      tib_xlim[2] = min(0.8, quantile(tib$cscore, probs = .99, na.rm = T))
    }

    #
    is_target_qval_filtered = tib$type == "target" & tib$detect
    target_cscore_min = min(tib$cscore[is_target_qval_filtered], na.rm = T)
    legend_text = sprintf("\n\nCscore: %.3f\n#pep: %d  %d%%", target_cscore_min, sum(is_target_qval_filtered), round(sum(is_target_qval_filtered) / length(is_target_qval_filtered) * 100))

    # plot
    p = suppressWarnings(ggplot(tib, aes(x = cscore, y = ..density.., fill = type)) +
                           geom_histogram(position = "identity", colour = "black", bins = 30, na.rm = T) +
                           scale_fill_manual(values = c("target" = "#01857188", "decoy" = "#aaaaaaff")) +
                           geom_vline(linetype = "dashed", color = "black", xintercept = target_cscore_min) +
                           annotate(geom = "text", label = legend_text, x = target_cscore_min, y = Inf, hjust = -0.1, vjust = 0.6, size = 2) +
                           scale_x_continuous(limits = tib_xlim, expand = c(0, 0)) +
                           labs(x = "", y = "", title = sid, fill = "") +
                           ggpubr::theme_classic2() +
                           theme(
                             legend.position = "none",
                             axis.line.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.text.y = element_blank(),
                             text = element_text(size = 9),
                             plot.title = element_text(size = 7, hjust = 0.5)
                           ))

    plotlist[[sid]] = p
    plotscores = c(plotscores, sum(is_target_qval_filtered))

    ## debug
    # hist(tib$cscore[tib$type=="target"], col="grey", freq=F)
    # hist(tib$cscore[tib$type=="decoy"], col="red", freq=F, add=T)
  }
  plotlist = plotlist[order(plotscores, decreasing = T, na.last = T)]
  rm(plotscores)

  result = list()
  # split all plots into sets of 3
  plot_list_by_row = split_array(plotlist, chunk_size = 3)
  for (l in plot_list_by_row) {
    # each set of 3 plots is a 'row' in output PDF, combine these plots into n columns
    result[[length(result) + 1]] = suppressWarnings(ggpubr::ggarrange(plotlist = l, ncol = 3, nrow = 1))
  }
  return(result)
}

