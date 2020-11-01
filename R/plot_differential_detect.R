
#' placeholder title
#' @param dataset dataset object
plot_differential_detect = function(dataset) {
  result = list()
  if(!"dd_proteins" %in% names(dataset) || nrow(dataset$dd_proteins) == 0) {
    return(result)
  }

  for(contr in dataset_contrasts(dataset)) {
    tib_contr = dataset$dd_proteins %>% filter(contrast == contr & !is.na(zscore_count_detect))
    if(nrow(tib_contr) == 0) next

    p_hist = ggplot(tib_contr, aes(zscore_count_detect)) +
      geom_histogram(bins=25, boundary = 0, colour = "white", fill="darkgrey", na.rm=T) +
      geom_vline(xintercept = c(-3, 3), colour = "red") +
      labs(x="differential detect z-score", y="number of proteins", colour = "", title = contr,
           subtitle = sprintf("#proteins tested: %d  #abs(zscore) >= 3: %d", nrow(tib_contr), sum(abs(tib_contr$zscore_count_detect)>=3, na.rm = T))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=8),
            plot.subtitle = element_text(hjust = 0.5, size=8),
            legend.position = "none")

    result[[contr]] = p_hist
  }

  return(result)
}
