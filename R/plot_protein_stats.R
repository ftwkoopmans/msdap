#' ### foldchange and standard error in each contrast, one plot for each protein
#' # ggplot is too slow for thousands of plots
#' # instead, use a vanilla R implementation (downside: more hardcoding and fiddling, upside: few seconds to plot 2000+ proteins)
#' #' @param dataset todo
#' #' @param pdf_file_path todo
#' #'
#' #' @export
#' plot_protein_stats = function(dataset, pdf_file_path) {
#'   tib_stats = stats_gather_results(dataset)
#'   if (nrow(tib_stats) == 0) {
#'     append_log("no stats have been computed so far, nothing to plot", type = "warning")
#'     return()
#'   }
#'
#'   append_log("plotting statistical results per protein...", type = "info")
#'
#'   tib_stats$notch_y_low = tib_stats$foldchange.log2 - (1.96 * tib_stats$standarderror)
#'   tib_stats$notch_y_high = tib_stats$foldchange.log2 + (1.96 * tib_stats$standarderror)
#'   # sort tibble to get most significant proteins first
#'   tib_stats = tib_stats %>% arrange(qvalue)
#'   uprot = unique(tib_stats$protein_id)
#'
#'   levels_methods = names(dataset$stats)
#'   levels_contrast = names(dataset$stats[[1]])
#'
#'   # setup pdf
#'   graphics.off() # bugfix
#'   Sys.sleep(1) # bugfix, sometimes opening/closing PDFs after one-another too quickly causes problems
#'   pdf(pdf_file_path, 2 + 3 * length(levels_methods) * length(levels_contrast), 6, pointsize = 11) # A4; 8.27, 11.69
#'   par(mfrow = c(1, length(levels_methods)), mar = c(9, 4, 4, 2), mgp = c(2, .75, 0)) # bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
#'   global_ylim = c(-2, 2) # default y-lim
#'   global_xlim = c(0.5, length(levels_contrast) + 0.5)
#'
#'   pdf_is_tiny = length(levels_methods) * length(levels_contrast) == 1
#'
#'
#'   for (prot in uprot) { # prot="IQSEC2"; prot=uprot[1]
#'     tib = tib_stats %>% filter(protein_id == prot)
#'
#'     # for the y-limits, use at least the global limit. if the current row is out-of-bounds, increase limit for current protein (exception / outlier)
#'     # note; same ylim between methods (eg; ebayes and msqrob)
#'     prot_ylim = c(
#'       floor(min(c(tib$notch_y_low, global_ylim), na.rm = T)),
#'       ceiling(max(c(tib$notch_y_high, global_ylim), na.rm = T))
#'     )
#'
#'     # iterate methods
#'     for (index_method in seq_along(levels_methods)) {
#'       name_method = levels_methods[index_method]
#'
#'
#'       plot(NA, type = "n", xlim = global_xlim, xaxt = "n", ylim = prot_ylim, las = 1, main = "", xlab = "", ylab = "log2 fold-change")
#'       if (index_method == 1) { # on first panel, print full protein name for reference
#'         mtext(paste(paste(unique(c(tib$protein_id[1], tib$gene_symbols_or_id[1])), collapse = " "), ifelse(pdf_is_tiny, "\n", ""), sub("^([^=]+) [A-Z]+=.*", "\\1 ...", tib$fasta_headers[1])),
#'           side = 3, line = 1.5, cex = ifelse(pdf_is_tiny, .5, .9), adj = 0
#'         ) # if we have a tiny PDF (eg; WT/KO, only eBayes), reduce title cex
#'       }
#'       mtext(name_method, side = 3, line = .5, cex = 1, adj = 0.5, font = 2) # method name
#'       # axis(side = 1, at=seq_along(levels_contrast), labels = levels_contrast, las=2)
#'       text(
#'         x = seq_along(levels_contrast), y = par()$usr[3] - 0.1 * (par()$usr[4] - par()$usr[3]), # https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/
#'         labels = levels_contrast, srt = 45, pos = 1, xpd = NA
#'       )
#'
#'       # iterate contrasts and plot
#'       for (index_contrast in seq_along(levels_contrast)) {
#'         name_contrast = levels_contrast[index_contrast]
#'         # subset data
#'         x = tib %>% filter(contrast_name == name_contrast & method_name == name_method)
#'         if (nrow(x) != 1) {
#'           next()
#'         }
#'
#'         # CI stick
#'         lines(
#'           x = c(index_contrast, index_contrast),
#'           y = c(x$notch_y_low, x$notch_y_high)
#'         )
#'         # CI notches
#'         lines(
#'           x = c(index_contrast - .1, index_contrast + .1),
#'           y = c(x$notch_y_low, x$notch_y_low)
#'         )
#'         lines(
#'           x = c(index_contrast - .1, index_contrast + .1),
#'           y = c(x$notch_y_high, x$notch_y_high)
#'         )
#'         # foldchange horizontal line
#'         lines(
#'           x = c(index_contrast - .25, index_contrast + .25),
#'           y = c(x$foldchange.log2, x$foldchange.log2), col = "blue", lwd = 2
#'         )
#'         # text
#'         text(x = index_contrast - .11, y = x$foldchange.log2, labels = sprintf("fc:%.2f", x$foldchange.log2), adj = c(1, -.5))
#'         text(x = index_contrast + .11, y = x$foldchange.log2, labels = sprintf("q:%.2e", x$qvalue), adj = c(0, -.5))
#'       }
#'     }
#'   } # end protein loop
#'
#'   dev.off()
#'   graphics.off()
#' }
#'
#'
#'
#' ### ggplot is just too slow for creating 1000's of plots
#' # plot_protein_stats = function(dataset) {
#' #   tib_stats = stats_gather_results(dataset)
#' #   if(nrow(tib_stats)==0) {
#' #     return(list())
#' #   }
#' #
#' #   ## prep data tibble
#' #   # convert to factors; this enforces order in downstream ggplot and together with drop=F guarantees each factor level is shown
#' #   tib_stats$contrast_name = factor(tib_stats$contrast_name, levels=names(dataset$stats[[1]]))
#' #   tib_stats$method_name = factor(tib_stats$method_name, levels=names(dataset$stats))
#' #   # prep plot data. do this on overall tibble to increase speed (as compared to doing within loop)
#' #   i = match(tib_stats$contrast_name, levels(tib_stats$contrast_name))
#' #   tib_stats$fc_left=i-.25
#' #   tib_stats$fc_right=i+.25
#' #   tib_stats$notch_left=i-.05
#' #   tib_stats$notch_right=i+.05
#' #   tib_stats$notch_y_low = tib_stats$foldchange.log2 - (1.96  * tib_stats$standarderror)
#' #   tib_stats$notch_y_high = tib_stats$foldchange.log2 + (1.96 * tib_stats$standarderror)
#' #   tib_stats$title = paste(tib_stats$protein_id, tib_stats$gene_symbols_or_id, sub("^([^=]+) [A-Z]+=.*", "\\1 ...", tib_stats$fasta_headers))
#' #
#' #   ## TODO: optionally, subset significant proteins (in either method_name)
#' #
#' #   # pdf("C:/temp/x.pdf", 8.27, 11.69, pointsize=11)
#' #   plotlist = list()
#' #
#' #   ## all ggplot functions and packages for wrapping over multiple pages are slow. tested; ggforce/ggplus/cowplot/ggpubr
#' #   ## instead, we'll manually define proteins for each page, subset tibble accordingly and create a single facet_grid
#' #   page_nrow = 4
#' #   page_ncol = 2
#' #   uprot = unique(tib_stats$protein_id)
#' #   nprot = length(uprot)
#' #   page_nprot = page_nrow * page_ncol
#' #   page_count = ceiling(nprot / page_nprot) # each page has ncol*nrow proteins
#' #   for(page_index in 1:page_count) { # page_index=104
#' #     ## proteins for current page  &  subset tibble
#' #     i = (page_index-1) * page_nprot
#' #     j = min(c(i+page_nprot-1, nprot))
#' #     prots = uprot[i:j]
#' #     # tib = tib_stats %>% filter(protein_id %in% prots)
#' #
#' #     ## ggplot each protein, store in list
#' #     l = list()
#' #     for(prot in prots) {
#' #       tib = tib_stats %>% filter(protein_id == prot)
#' #       # l[[prot]]
#' #       plotlist[[length(plotlist) + 1]] = ggplot(tib, aes(x=contrast_name, y=foldchange.log2), drop=F) +
#' #         geom_segment(aes(x = contrast_name, y = notch_y_low, xend = contrast_name, yend = notch_y_high)) +
#' #         geom_segment(aes(x = notch_left, y = notch_y_low, xend = notch_right, yend = notch_y_low)) +
#' #         geom_segment(aes(x = notch_left, y = notch_y_high, xend = notch_right, yend = notch_y_high)) +
#' #         geom_segment(aes(x = fc_left, y = foldchange.log2, xend = fc_right, yend = foldchange.log2), colour="blue", size=1.2) +
#' #         geom_text(aes(label=sprintf("fc:%.2f", foldchange.log2)), hjust=1.25, vjust=-1) +
#' #         geom_text(aes(label=sprintf("p:%.1e", qvalue)), hjust=-0.25, vjust=-1) +
#' #         expand_limits(y=c(-2,2)) +
#' #         facet_grid(~method_name, drop = F) +
#' #         labs(x="", y="log2 fold-change", title=tib$title[1]) +
#' #         theme_bw() +
#' #         theme(line = element_blank(), title = element_text(size=6), axis.title.y = element_text(size=8))
#' #     }
#' #     ## combine ggplot of each protein using ggpubr::ggarrange
#' #     # plotlist[[length(plotlist) + 1]] = ggarrange(plotlist = l, nrow=page_nrow, ncol=page_ncol)
#' #     # print(ggarrange(plotlist = l, nrow=page_nrow, ncol=page_ncol))
#' #   }
#' #   # dev.off()
#' #
#' #   return(plotlist)
#' # }
#'
#' # t1 <- Sys.time()
#' # l = plot_protein_stats(dataset)
#' # t2 <- Sys.time()
#' # print(difftime(t2, t1, units = "secs"))
#' #
#' # t1 <- Sys.time()
#' # pdf("C:/temp/x.pdf", 8.27, 11.69, pointsize=11)
#' # l
#' # dev.off()
#' # t2 <- Sys.time()
#' # print(difftime(t2, t1, units = "secs"))
