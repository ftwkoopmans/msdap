
#' placeholder title
#' @param dataset dataset object
plot_differential_detect = function(dataset) {

  tib = dataset$dd_proteins %>%
    left_join(dataset$proteins %>% select(protein_id, gene_symbols_or_id), by="protein_id") %>%
    filter(!is.na(diff_detect_zscore)) %>%
    arrange(desc(abs(diff_detect_zscore)))

  result = list()

  for(contr in unique(tib$contrast)) {
    tib_contr = tib %>% filter(contrast == contr)
    # x = tib_contr$diff_detect_zscore[tib_contr$diff_detect_zscore_candidate]

    p_hist = ggplot(tib_contr, aes(diff_detect_zscore)) +
      geom_histogram(bins=25, boundary = 0, colour = "white", fill="darkgrey", na.rm=T) +
      geom_vline(xintercept = c(-2, 2), colour = "red") +
      # geom_vline(xintercept = c(max(x[x<0]), min(x[x>1])), colour = "red") +
      labs(x="differential detect z-score", y="number of proteins", colour = "", title = contr, subtitle = sprintf("#proteins:%d #candidates:%d", nrow(tib_contr), sum(tib_contr$diff_detect_zscore_candidate))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=10),
            plot.subtitle = element_text(hjust = 0.5, size=8),
            legend.position = "none")

    result[[contr]] = p_hist

    # g = tib_contr %>% filter(diff_detect_zscore_candidate) %>% pull(gene_symbols_or_id)
    # g = sapply(g, function(x) paste(stringr::str_trunc(x, 12, "right")) )
    # result[[contr]] = list(plot = p_hist, genes = g)

  }

  return(result)

  # tib_input = dataset$dd_proteins %>% filter(contrast == dataset$dd_proteins$contrast[1])
  # x = tib_input$diff_detect_zscore[dataset$dd_proteins$diff_detect_zscore_candidate]
  # hist(tib_input$diff_detect_zscore, breaks=25)
  # abline(v=c(max(x[x<0]), min(x[x>1])), col=2)
}
