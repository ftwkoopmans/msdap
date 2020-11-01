
#' placeholder title
#' @param stats_de todo
#' @param log2foldchange_threshold todo
#' @param qvalue_threshold todo
#' @param mtitle todo
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggrepel geom_text_repel
plot_volcano = function(stats_de, log2foldchange_threshold = NA, qvalue_threshold=NA, mtitle = "") {

  ########### format input data
  stats_de$pvalue[!is.finite(stats_de$pvalue)] = NA
  stats_de$qvalue[!is.finite(stats_de$qvalue)] = NA

  # optionally, if the user wants to re-define the qvalue cutoff for this plot, do that first. Next, optionally take all significant hits and add filtering by foldchange
  if(!is.na(qvalue_threshold)) {
    stats_de = stats_de %>% mutate(signif = is.finite(qvalue) & qvalue <= qvalue_threshold)
  }
  if(!is.na(log2foldchange_threshold)) {
    log2foldchange_threshold = abs(log2foldchange_threshold)
    stats_de = stats_de %>% mutate(signif = signif & abs(foldchange.log2) >= log2foldchange_threshold)
  }

  ## pretty-print labels
  stats_de$label = stats_de$protein_id
  # use gene symbols if available
  if("gene_symbols_or_id" %in% colnames(stats_de)) {
    stats_de$label = stats_de$gene_symbols_or_id
  }
  # strip ambiguous IDs
  stats_de$label = gsub(";.*", "", stats_de$label)
  # reduce the string length of very long labels
  rows = nchar(stats_de$label) > 10
  stats_de$label[rows] = paste0(substr(stats_de$label[rows], 1, 7), "...")
  # mark ambiguous IDs with an asterix
  rows = grepl(";", stats_de$label, fixed=T)
  stats_de$label[rows] = paste0(stats_de$label[rows], "*")



  ########### create volcano plots

  volcano_plotlist = list()
  for (algo_name in unique(stats_de$algo_de)) { #algo_name = "ebayes"
    # prepare data for current contrast
    tib = stats_de %>%
      filter(algo_de == algo_name) %>%
      drop_na(foldchange.log2, pvalue, qvalue) %>%
      # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
      mutate(minlog10qval = minlog10(qvalue),
             foldchange.log2_abs = abs(foldchange.log2)) %>%
      # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
      arrange(desc(pvalue)) %>%
      # reduce tibble size
      select(label, foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

    # which proteins should get a text label? all significant is not an option, some tests have hundreds
    # if |signif| < x, x, otherwise 10?  or flat amount?
    tib$flag_plot_label = rep(c(F,T), c(nrow(tib)-25, 25)) # plot top25, assuming data was sorted by pvalue descending

    # classify up/down regulated
    tib$updown = ifelse(tib$foldchange.log2 < 0, "down", "up")
    tib$updown[tib$signif != TRUE] = "unchanged"

    ### find outliers, values that are so far away that they may skew the plot's appearance, and classify data points accordingly
    xmax_nooutlier = c(-1,1) * max(tib$foldchange.log2_abs, na.rm = T)
    ymax_nooutlier = c(0, max(tib$minlog10qval, na.rm=T))
    xmax = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(.005, .995), na.rm = T)))
    ymax = c(0, quantile(tib$minlog10qval, probs = .995, na.rm = T))

    tib$outlier_x_low = tib$foldchange.log2 < xmax[1]
    tib$outlier_x_high = tib$foldchange.log2 > xmax[2]
    tib$outlier_y = tib$minlog10qval > ymax[2]

    tib$x_outlier = tib$foldchange.log2
    tib$y_outlier = tib$minlog10qval
    tib$x_outlier[tib$outlier_x_low] = xmax[1]
    tib$x_outlier[tib$outlier_x_high] = xmax[2]
    tib$y_outlier[tib$outlier_y] = ymax[2]

    tib$updown_outlier = tib$updown
    tib$updown_outlier[tib$updown == "unchanged" & (tib$outlier_x_low | tib$outlier_x_high)] = "unchanged_outlier"
    tib$updown_outlier[tib$updown == "down" & (tib$outlier_x_low | tib$outlier_y)] = "down_outlier"
    tib$updown_outlier[tib$updown == "up" & (tib$outlier_x_high | tib$outlier_y)] = "up_outlier"

    ### construct facets
    plottype_labels = c(asis = "data as−is, no labels",
                        asis_lab = "data as−is, label 25 best qvalue",
                        lim = "limited x− and y−axis, no labels",
                        lim_lab = "limited x− and y−axis, label 25 best qvalue")

    tib_facets = bind_rows(tib %>% select(label, x=foldchange.log2, y=minlog10qval, pch=updown, flag_plot_label) %>% mutate(flag_plot_label=FALSE) %>% add_column(plottype = "asis"),
                           tib %>% select(label, x=foldchange.log2, y=minlog10qval, pch=updown, flag_plot_label)                                   %>% add_column(plottype = "asis_lab"),
                           tib %>% select(label, x=x_outlier, y=y_outlier, pch=updown_outlier, flag_plot_label) %>% mutate(flag_plot_label=FALSE)  %>% add_column(plottype = "lim"),
                           tib %>% select(label, x=x_outlier, y=y_outlier, pch=updown_outlier, flag_plot_label)                                    %>% add_column(plottype = "lim_lab"))
    tib_facets$pch = factor(tib_facets$pch, levels = c("down", "unchanged", "up", "down_outlier", "unchanged_outlier", "up_outlier"))

    # some mock data to enforce symmetric x-axis and expand the y-limit a bit to make room for the labels (this is a workaround because we cannot hardcode separate y-axis limits per facet)
    blank_data = bind_rows(tibble(x=xmax_nooutlier, y=ymax_nooutlier[2] * 1.2, plottype = "asis"),
                           tibble(x=xmax_nooutlier, y=ymax_nooutlier[2] * 1.2, plottype = "asis_lab"),
                           tibble(x=xmax, y=ymax[2] * 1.2, plottype = "lim"),
                           tibble(x=xmax, y=ymax[2] * 1.2, plottype = "lim_lab") ) %>%
      add_column(pch="unchanged", label="")

    ### volcano plot
    p = ggplot(tib_facets, aes(x, y, colour = pch, fill = pch, shape = pch, label = label)) +
      geom_point(na.rm=T) +
      geom_blank(data = blank_data) +
      scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        values = c(up = "#d55e00aa", down = "#56b4e9aa", up_outlier = "#d55e00aa", down_outlier = "#56b4e9aa", unchanged = "#22222299", unchanged_outlier = "#22222299"),
        labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
                   unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
        drop = F,
        guide = guide_legend(byrow=T, override.aes = list(alpha = 1))
      ) +
      scale_shape_manual(values = c(up = 24, down = 25, up_outlier = 2, down_outlier = 6, unchanged = 19, unchanged_outlier = 1),
                         labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
                                    unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
                         drop = F,
                         guide = guide_legend(byrow=T, override.aes = list(alpha = 1))
      ) +
      ggrepel::geom_text_repel(alpha=1, data = tib_facets %>% filter(flag_plot_label == TRUE), segment.alpha = 0.3, min.segment.length = unit(0.25, 'lines'), vjust = 0.6, show.legend = FALSE, size = 2) +
      facet_wrap(~plottype, nrow = 2, ncol = 2, scales = "free", labeller = labeller(plottype=plottype_labels)) +
      labs(x = "log2 fold-change", y = "-log10 FDR adjusted p-value", title = paste(algo_name, "@", mtitle)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=8),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size=8))

    if(!is.na(log2foldchange_threshold) & log2foldchange_threshold != 0) {
      p = p + geom_vline(xintercept = c(-1,1) * log2foldchange_threshold, colour = "darkgrey", linetype = "dashed")
    }
    if(!is.na(qvalue_threshold) & qvalue_threshold != 0) {
      p = p + geom_hline(yintercept = -log10(qvalue_threshold), colour = "darkgrey", linetype = "dashed")
    }

    volcano_plotlist[[algo_name]] = p
  }
  return(volcano_plotlist)
}
#### v1 below
# plot_volcano = function(stats_de, log2foldchange_threshold = NA, qvalue_threshold=NA, mtitle = "") {
#   stats_de$pvalue[!is.finite(stats_de$pvalue)] = NA
#   stats_de$qvalue[!is.finite(stats_de$qvalue)] = NA
#
#   # optionally, if the user wants to re-define the qvalue cutoff for this plot, do that first. Next, optionally take all significant hits and add filtering by foldchange
#   if(!is.na(qvalue_threshold)) {
#     stats_de = stats_de %>% mutate(signif = is.finite(qvalue) & qvalue <= qvalue_threshold)
#   }
#   if(!is.na(log2foldchange_threshold)) {
#     log2foldchange_threshold = abs(log2foldchange_threshold)
#     stats_de = stats_de %>% mutate(signif = signif & abs(foldchange.log2) >= log2foldchange_threshold)
#   }
#
#   volcano_plotlist = list()
#   for (algo_name in unique(stats_de$algo_de)) { #algo_name = "ebayes"
#     # order by p-value, draw most significant last
#     tib = stats_de %>% filter(algo_de == algo_name) %>% drop_na(foldchange.log2, pvalue, qvalue) %>% arrange(desc(pvalue))
#     tib$minlog10qval = minlog10(tib$qvalue)
#     # some stat methods, such as msempire, could give a zero for pval or qval. set -log10 to best non-finite value we got
#     # tib$minlog10qval[is.finite(tib$qvalue) & tib$qvalue == 0] = max(tib$minlog10qval[is.finite(tib$minlog10qval)])
#     # tib = tib %>% filter(is.finite(minlog10qval))
#
#     tib$updown = factor(ifelse(tib$foldchange.log2 < 0, "down", "up"), levels = c("up", "down", "up_outlier", "down_outlier", "unchanged", "unchanged_outlier"))
#     tib$updown[tib$signif != TRUE] = "unchanged"
#
#     # which proteins should get a text label? all significant is not an option, some tests have hundreds
#     # if |signif| < x, x, otherwise 10?  or flat amount?
#     tib$flag_plot_label = rep(c(F,T), c(nrow(tib)-25, 25)) # plot top25
#
#     # find outliers
#     ymax = c(0, ceiling(quantile(tib$minlog10qval, probs = .98, na.rm = T)))
#     xmax = c(-1,1) * max(ceiling(quantile(tib$foldchange.log2, probs = c(.02, .98), na.rm = T)))
#     tib$outlier_x_low = tib$foldchange.log2 < xmax[1]
#     tib$outlier_x_high = tib$foldchange.log2 > xmax[2]
#     tib$outlier_y = tib$minlog10qval > ymax[2]
#
#     tib$x_outlier = tib$foldchange.log2
#     tib$y_outlier = tib$minlog10qval
#     tib$x_outlier[tib$outlier_x_low] = xmax[1]
#     tib$x_outlier[tib$outlier_x_high] = xmax[2]
#     tib$y_outlier[tib$outlier_y] = ymax[2]
#     # tib$updown[tib$signif != 0 & (tib$outlier_x_low | tib$outlier_x_high)] = ""
#     tib$updown_outlier = tib$updown
#     tib$updown_outlier[tib$updown == "unchanged" & (tib$outlier_x_low | tib$outlier_x_high)] = "unchanged_outlier"
#     tib$updown_outlier[tib$updown == "down" & (tib$outlier_x_low | tib$outlier_y)] = "down_outlier"
#     tib$updown_outlier[tib$updown == "up" & (tib$outlier_x_high | tib$outlier_y)] = "up_outlier"
#
#     #
#     volcano_baseplot = function(tib, log2foldchange_threshold=NA, qvalue_threshold=NA) {
#       p = ggplot(tib, aes(x, y, colour = updown, fill = updown, shape = updown, label = label)) +
#         geom_point(na.rm=T) +
#         scale_discrete_manual(
#           aesthetics = c("colour", "fill"),
#           values = c(up = "#d55e00bb", down = "#56b4e9bb", up_outlier = "#d55e00bb", down_outlier = "#56b4e9bb", unchanged = "#000000bb", unchanged_outlier = "#000000bb"),
#           labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
#                      unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
#           drop = F
#         ) +
#         scale_shape_manual(values = c(up = 24, down = 25, up_outlier = 2, down_outlier = 6, unchanged = 19, unchanged_outlier = 1),
#                            labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
#                                       unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
#                            drop = F) +
#         labs(x = "log2 fold-change", y = "-log10 FDR adjusted p-value", title = paste(algo_name, "@", mtitle)) +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5, size=8),
#               plot.subtitle = element_text(hjust = 0.5, size=7),
#               legend.title = element_blank()) # legend.position = "none",
#
#       if(!is.na(log2foldchange_threshold)) {
#         p = p + geom_vline(xintercept = c(-1,1) * log2foldchange_threshold, colour = "darkgrey", linetype = "dashed")
#       }
#       if(!is.na(qvalue_threshold)) {
#         p = p + geom_hline(yintercept = -log10(qvalue_threshold), colour = "darkgrey", linetype = "dashed")
#       }
#
#       return(p)
#     }
#
#     tib$label = tib$protein_id
#     if("gene_symbols_or_id" %in% colnames(tib)) {
#       tib$label = tib$gene_symbols_or_id
#     }
#     tib$label = gsub(";.*", "", tib$label)
#     tib$label = ifelse(nchar(tib$label) <= 10, tib$label, paste0(substr(tib$label, 1, 7), "..."))
#
#     tib_plot_normal = tib %>% rename(x=foldchange.log2, y=minlog10qval)
#     tib_plot_outlier = tib %>% select(-updown) %>% rename(x=x_outlier, y=y_outlier, updown=updown_outlier)
#
#     p_volcano = volcano_baseplot(tib_plot_normal, log2foldchange_threshold, qvalue_threshold) +
#       ggplot2::xlim(c(-1,1) * max(abs(tib_plot_normal$x), na.rm=T)) +
#       ggplot2::ylim(c(0, max(tib_plot_normal$y, na.rm=T) * 1.2))
#
#     p_volcano_outliers_limited = volcano_baseplot(tib_plot_outlier, log2foldchange_threshold, qvalue_threshold) +
#       ggplot2::xlim(xmax) +
#       ggplot2::ylim(c(0, ymax[2] * 1.2))
#
#     ## add text labels
#     p_volcano__labeled = p_volcano +
#       ggrepel::geom_text_repel(data = tib_plot_normal %>% mutate(label = ifelse(flag_plot_label, label, "")), segment.alpha = 0.4, min.segment.length = unit(0.25, 'lines'), vjust = 0.6, show.legend = FALSE, size = 2)
#
#     p_volcano_outliers_limited__labeled = p_volcano_outliers_limited +
#       ggrepel::geom_text_repel(data = tib_plot_outlier %>% mutate(label = ifelse(flag_plot_label, label, "")), segment.alpha = 0.4, min.segment.length = unit(0.25, 'lines'), vjust = 0.6, show.legend = FALSE, size = 2)
#
#     volcano_plotlist[[algo_name]] = ggpubr::ggarrange(plotlist = list(p_volcano + labs(subtitle = "data as-is, no labels"),
#                                                                       p_volcano__labeled + labs(subtitle = "data as-is, label 25 best qvalue"),
#                                                                       p_volcano_outliers_limited + labs(subtitle = "limited x- and y-axis, no labels"),
#                                                                       p_volcano_outliers_limited__labeled + labs(subtitle = "limited x- and y-axis, label 25 best qvalue")
#                                                                       ), ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
#   }
#
#   return(volcano_plotlist)
# }



#' placeholder title
#' @param tib pre-define color_code column
#' @param mtitle todo
plot_foldchanges = function(tib, mtitle="") {
  ggplot(tib, aes(foldchange.log2, colour = color_code)) +
    geom_vline(xintercept=0) +
    # same density function as default in base R's stats::density()
    geom_line(stat = "density", adjust=1, bw = "SJ", na.rm = T) +
    facet_wrap( ~ algo_de, nrow = ceiling(sqrt(n_distinct(tib$color_code))), scales = "free_y") +
    coord_cartesian(xlim = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(0.005, 0.995), na.rm=T)))) +
    labs(x="log2 foldchange", colour = "", title=mtitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=8),
          legend.position = "none") # bottom

  ### alternatively, combine into a single plot. Since some methods (eg; msqrob) can hugely inflate the amount of proteins with exactly zero foldchange this can distort the binning and y-axis scale. So above plot with separate facets (and ditto free scales) is preferred
  # ggplot(tib, aes(x = foldchange.log2, colour=color_code)) +
  #   geom_vline(xintercept=0) +
  #   stat_density(adjust=3, bw = "SJ", geom="line", na.rm=T, trim=F) + # aes(y = ..count..) seems scuffed, way to many
  #   coord_cartesian(xlim = c(-1,1) * mean(abs(quantile(tib$foldchange.log2, probs = c(0.005, 0.995), na.rm=T)))) +
  #   labs(x="log2 foldchange", colour = "", title=mtitle) +
  #   theme_bw() +
  #   theme(legend.position = "bottom")
}



#' placeholder title
#' @param tib todo
#' @param mtitle todo
plot_pvalue_histogram = function(tib, mtitle="") {
  # debug_tib_pvalue_hist <<- tib
  binwidth = 0.05
  tib_summ = tib %>% filter(is.finite(pvalue)) %>% group_by(algo_de) %>% summarise(yintercept_line = binwidth * n())

  ggplot(tib, aes(pvalue)) +
    geom_histogram(binwidth = binwidth, boundary = 0, colour = "white", fill="darkgrey", na.rm=T) +
    geom_hline(data = tib_summ, aes(yintercept = yintercept_line), colour = "darkblue") +
    scale_x_continuous(breaks = (0:10)/10, minor_breaks = (0:5)/5) +
    facet_wrap( ~ algo_de, nrow = 2, scales = "fixed") +
    labs(x="p-value", y="number of proteins", colour = "", title=mtitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=8),
          legend.position = "none")
}
