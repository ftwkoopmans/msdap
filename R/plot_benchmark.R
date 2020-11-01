
#' placeholder title
#' prior to this; species %in% c("YEAST", "HUMAN")
#' add extra column predictor; true for spike-in, false for rest
#' @param tib todo
#' @param mtitle todo
#' @param universe todo
#' @param plot_coords todo
#'
#' @importFrom pROC plot.roc coords
#' @importFrom gtools mixedsort
#' @importFrom colorspace hcl_palettes
plot_roc = function(tib, mtitle="", universe="all", plot_coords = FALSE) {
  if(!(universe %in% c("all", "signif", "signif_any"))) {
    append_log('universe parameter must be any of; "all", "signif", "signif_any"', type = "error")
  }

  oldpar = par(mfrow = c(2, 2))
  ualg = mixedsort(unique(tib %>% pull(algo_de)))
  # enforce sorting; if only plotting DE algorithms, first show protein-level models then peptide-level models. otherwise, just string sort
  if(all(ualg %in% c("ebayes", "msqrobsum", "msempire", "msqrob"))) {
    ualg = ualg[order(match(ualg, c("ebayes", "msqrobsum", "msempire", "msqrob")))]
  } else {
    ualg = mixedsort(ualg)
  }
  # for padding legend downstream
  ualg_max_nchar = max(nchar(ualg))

  # clrs = paste0(brewer.pal(n = 12, name = 'Paired'), "BB")
  # clrs = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
  clrs = colorspace::hcl_palettes(n = 9, palette = 'Dark3')

  qval_cutoff = 0.01
  min_fc=1
  # for (min_fc in c(1, 1.1)) { # min_fc=1.1
  for (plot_coords in c(F, T)) { # qval_cutoff = 0.05
    # for (qval_cutoff in c(0.01, 0.05)) { # qval_cutoff = 0.05
    for (zoom_plot in c(F, T)) { # qval_cutoff = 0.05
      # the universe for the ROC plot is all proteins that match filter rules in any method/analysis
      if(universe == "all") {
        tib_valid_protein_id = tib %>%
          filter(abs(foldchange.log2) >= log2(min_fc)) %>%
          select(protein_id, predictor) %>%
          distinct(protein_id, .keep_all = T)
      }
      if(universe == "signif_any") {
        tib_valid_protein_id = tib %>%
          filter(qvalue <= qval_cutoff & abs(foldchange.log2) >= log2(min_fc)) %>%
          select(protein_id, predictor) %>%
          distinct(protein_id, .keep_all = T)
      }

      # alg_lty = 1 + ((seq_along(ualg)-1) %% 4)
      # alg_lty = ceiling(seq_along(ualg) / 12)

      alg_lty = rep(1, length(ualg))
      # if(length(ualg) > 5) {
        # alg_lty = c(1, 2, 5)[1 + ((seq_along(ualg)-1) %% 3)]
      # }
      alg_clr = clrs[1 + ((seq_along(ualg)-1) %% 12)]
      lgnd = NULL
      roc_plot_drawn = FALSE
      for (index_alg in seq_along(ualg)) { # index_alg = 1
        alg = ualg[index_alg]
        if(universe == "signif") {
          tib_valid_protein_id = tib %>%
            filter(algo_de == alg & qvalue <= qval_cutoff & abs(foldchange.log2) >= log2(min_fc)) %>%
            select(protein_id, predictor) %>%
            distinct(protein_id, .keep_all = T)
        }

        x = tib %>%
          filter(algo_de == alg & protein_id %in% tib_valid_protein_id$protein_id & abs(foldchange.log2) >= log2(min_fc)) %>%
          select(protein_id, pvalue, predictor, qvalue)
        # suppose that some of the proteins of interest are not in results from current method, we add them with a poor score
        if (any(!tib_valid_protein_id$protein_id %in% x$protein_id)) {
          x = rbind(x,
                    tib_valid_protein_id %>% filter(!protein_id %in% x$protein_id) %>% add_column(pvalue = 1, qvalue = 1, .after = "protein_id"))
        }

        roc_predictor_score = x$pvalue
        roc_classification = x$predictor
        if(nrow(x) == 0 || length(unique(roc_classification)) != 2) {
          next
        }

        # @ pROC documentation
        # "levels = the value of the response for controls and cases respectively"
        # "direction = in  which  direction  to  make  the  comparison? “>”: if the predictor values for the control group are higher than the values of thecase group (controls > t >= cases"
        # -->> since we ROC the p-values, control values are higher than case values
        roc_obj = pROC::plot.roc(roc_classification, roc_predictor_score, levels=c("background", "foreground"), direction=">", partial.auc = c(1, 0.9),
                                 xlim=c(1, ifelse(zoom_plot, .8, 0)), ylim=c(ifelse(zoom_plot, .5, 0), 1),
                              col = paste0(alg_clr[index_alg], "BB"), lty = alg_lty[index_alg], lwd=1,
                              main = sprintf("%s universe:%s FC>=%s qval<=%s", mtitle, universe, min_fc, qval_cutoff), cex.main = .6, add=roc_plot_drawn)

        coord_best = pROC::coords(roc_obj, "best", ret="all", transpose = FALSE, best.method="youden")
        if(plot_coords) {
          points(coord_best$specificity, coord_best$sensitivity, col=alg_clr[index_alg], pch="|", cex=1.5) # pch = 3 = +
          coord_signif = coords(roc_obj, 0.01, transpose = FALSE)
          points(coord_signif$specificity, coord_signif$sensitivity, col=alg_clr[index_alg], pch="\\", cex=1.5) # pch = 3 = x
        }

        lgnd = c(lgnd, sprintf("%s #:%d qval<=%s:%d TP:%d FP:%d pAUC:%.1f%% bestT:%.4g",
                               paste0(alg, paste(rep("  ", ualg_max_nchar - nchar(alg)), collapse = "")),  # simple string padding by double space
                               nrow(tib_valid_protein_id), qval_cutoff, sum(x$qvalue <= qval_cutoff), sum(x$qvalue <= qval_cutoff & x$predictor == "foreground"), sum(x$qvalue <= qval_cutoff & x$predictor == "background"), as.numeric(roc_obj$auc) * 100,
                               coord_best$threshold))
        roc_plot_drawn = T
      }

      if(roc_plot_drawn & !zoom_plot) {
        abline(v = .9, col = "grey", lty = 3)
        legend("bottomright", legend = lgnd, col = alg_clr, lty = alg_lty, lwd = 1.5, cex = .55, bty = "n")
      }
    }
  }
  par(oldpar)
}



#' placeholder title
#' pre-define color_code and minlog10pval
#' @param tib todo
#' @param mtitle todo
#' @param color_discrete todo
#'
#' @importFrom viridis scale_fill_viridis
plot_benchmark_volcano = function(tib, mtitle="", color_discrete = TRUE) {
  xlim_symmetric = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(0.001, 0.999), na.rm=T)))
  # xlim_symmetric = c(-1,1) * max(tib$foldchange.log2, na.rm=T)

  if(!color_discrete) {
    tib = tib %>% arrange(color_code)
  }

  p = ggplot(tib, aes(foldchange.log2, minlog10pval, colour = color_code, fill = color_code))
  if("shape" %in% colnames(tib)) {
    p = p + geom_point(alpha = 0.5, aes(shape=shape)) #+ scale_shape(solid = FALSE)
  } else {
    p = p + geom_point(alpha = 0.5)
  }
  p = p +
    facet_wrap( ~ algo_de, ncol = 2, scales = "free_y") +
    xlim(xlim_symmetric) +
    labs(x = "log2 foldchange (far extremes not shown)", y = "-log10 FDR adjusted p-value", title = mtitle) +
    theme_bw() +
    theme(legend.position = "bottom")
  if(color_discrete) {
    p = p + scale_color_brewer(palette="Dark2")
  } else {
    # p = p + viridis::scale_fill_viridis(option = "E", aesthetics = c("colour", "fill"))
    p = p + scale_fill_distiller(palette = "Spectral", aesthetics = c("colour", "fill"))
  }

  return(p)
}



#' placeholder title
#' @param tib todo
#' @param mtitle todo
#'
#' @importFrom ggpubr ggarrange
plot_true_false_positive_counts = function(tib, mtitle = "") {
  plotlist = list()
  for (min_fc in c(1, 1.1)) { # min_fc=1.1
    lbl = sprintf("%s FC>=%s", mtitle, min_fc)
    # bins = seq(0.01, 0.1, by=0.01)
    bins = seq(0.001, 0.05, by=0.001)

    x = tib %>%
      select(protein_id, qvalue, predictor, algo_de, foldchange.log2) %>%
      filter(qvalue <= max(bins) & abs(foldchange.log2) >= log2(min_fc))

    x$bin = 1
    for(i in 1:length(bins)) { #i=4
      rows = x$qvalue > bins[i] # & x$qvalue > bins[i+1]
      x$bin[rows] = i
    }

    x2 = tibble()
    for(alg in unique(x$algo_de)) {
      for(pred in unique(x$predictor)) {
        for(b in seq_along(bins)) {
          tib_subset = x %>% filter(algo_de==alg & predictor == pred & bin <= b)
          x2 = bind_rows(x2, tibble(algo_de=alg, predictor=pred, bin=b, cs = nrow(tib_subset)))
        }
      }
    }
    # dplyr summary would skip the bins that are the same
    # x2 = x %>% count(algo_de, predictor, bin) %>% group_by(algo_de, predictor) %>% mutate(cs = cumsum(n)) %>% ungroup()


    x2$predictor = factor(x2$predictor, levels=c("foreground", "background"))
    x2$bin_qval = bins[x2$bin]

    p = ggplot(x2, aes(x=bin_qval, y=cs, fill = algo_de, colour = algo_de)) +
      # geom_step(size=1.5) + facet_wrap(~predictor) +
      geom_step(aes(linetype=predictor)) + #size=1.5
      # scale_x_continuous(breaks = (0:10 + .5)/100, labels = (0:10)/100, minor_breaks = NULL) +
      labs(title=lbl, x = "Q-value threshold (bin)", y="Cumulative protein count") +
      theme_bw()
    plotlist[[length(plotlist) + 1]] = p

    # FPR = fraction of cumsum; sum(background) / sum(foreground+background)
    x_fdr = x2 %>%
      group_by(algo_de, bin, bin_qval) %>%
      add_tally(name = "n_total", wt = cs) %>%
      filter(predictor == "background") %>%
      summarise(fpr = cs / n_total * 100)

    p = ggplot(x_fdr, aes(x=bin_qval, y=fpr, fill = algo_de, colour = algo_de)) +
      geom_step() + # size=1.5
      # scale_x_continuous(breaks = (0:10 + .5)/100, labels = (0:10)/100, minor_breaks = NULL) +
      labs(title=lbl, x = "Q-value threshold (bin)", y="False Positive Rate for proteins at some confidence threshold (%)") +
      theme_bw()
    plotlist[[length(plotlist) + 1]] = p
  }

  ggpubr::ggarrange(plotlist = plotlist, ncol=2, nrow=2)
}
