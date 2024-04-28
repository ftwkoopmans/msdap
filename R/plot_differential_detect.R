
#' Plot differential detection results as a histogram
#'
#' @param dataset dataset object
#' @param zscore_threshold cutoff used in the plot (for absolute values)
#' @returns list of ggplot objects with 1 plot per contrast. If no differential detection data is available, returns an empty list
#' @export
plot_differential_detect = function(dataset, zscore_threshold = 4) {
  result = list()
  if(!"dd_proteins" %in% names(dataset) || nrow(dataset$dd_proteins) == 0) {
    return(result)
  }

  for(contr in dataset_contrasts(dataset)) {
    tib_contr = dataset$dd_proteins %>% filter(contrast == contr & is.finite(zscore))
    if(nrow(tib_contr) == 0) {
      next
    }

    tib_contr = tib_contr %>% filter(type %in% c("detect", "quant"))
    tib_contr_detect = tib_contr %>% filter(type == "detect")
    tib_contr_quant = tib_contr %>% filter(type == "quant")
    lbl_detect = sprintf("only detected peptides; #proteins tested: %d  #abs(zscore) >= %s: %d", nrow(tib_contr_detect), as.character(zscore_threshold), sum(abs(tib_contr_detect$zscore) >= zscore_threshold))
    lbl_quant = sprintf("all quantified peptides; #proteins tested: %d  #abs(zscore) >= %s: %d", nrow(tib_contr_quant), as.character(zscore_threshold), sum(abs(tib_contr_quant$zscore) >= zscore_threshold))

    tib_contr = tib_contr %>%
      arrange(type) %>%
      mutate(type_label = ifelse(type == "detect", lbl_detect, lbl_quant),
             type_label = factor(type_label, levels = unique(type_label)))

    p_hist = ggplot(tib_contr, aes(zscore)) +
      geom_histogram(bins=25, boundary = 0, colour = "black", fill="lightgrey", na.rm=T) +
      geom_vline(xintercept = c(-zscore_threshold, zscore_threshold), colour = "red") +
      facet_wrap(.~type_label, ncol = 1, scales = "free") +
      labs(x="Differential z-score for observed peptides", y="Protein count", colour = "", title = contr) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=9),
            plot.subtitle = element_text(hjust = 0.5, size=9),
            legend.position = "none")

    result[[contr]] = p_hist
  }

  return(result)
}
