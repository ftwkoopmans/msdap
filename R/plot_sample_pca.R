
#' placeholder title
#' for reference, standard pca (does not cope with missing values); summary(stats::prcomp(t(matrix_sample_intensities), center = T, scale. = F))
#' @param matrix_sample_intensities todo
#' @param samples todo
#' @param samples_colors todo
#' @param label_samples_by_shortname todo
#' @param pch_as_exclude todo
#'
#' @importFrom pcaMethods pca
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr ggarrange
#' @importFrom gtools mixedsort
plot_sample_pca = function(matrix_sample_intensities, samples, samples_colors, label_samples_by_shortname = TRUE, pch_as_exclude = TRUE) {
  PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
  mat_pca = PPCA@scores # rownames = sample_id
  pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores)))

  props = grep("^(sample_id$|sample_index$|shortname$|contrast:)", colnames(samples_colors), ignore.case = T, value = T, invert = T)

  # how should we label the samples? if not by shortname, use the index (and fallback to sample_id if index not available)
  if(length(label_samples_by_shortname) != 1 || !is.logical(label_samples_by_shortname) || !is.finite(label_samples_by_shortname)) {
    label_samples_by_shortname = nrow(samples) < 20
  }
  prop_sample_labels = "shortname"
  if(!label_samples_by_shortname) {
    prop_sample_labels = intersect(c("sample_index", "sample_id"), colnames(samples))[1]
  }

  pcaplots = list()
  for (prop in props) { #prop=props[1]
    # don't plot 'exclude' property if there are no excluded samples
    if(prop == "exclude" && !any(samples %>% filter(sample_id %in% rownames(mat_pca)) %>% pull(exclude))) {
      next
    }

    plotlist = list()
    for (dims in list(1:2, c(1, 3), 2:3)) { # dims=1:2
      # extract data from PCA object  &  join with sample metadata to get color-coding and label
      tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>%
        left_join(samples %>% select(sample_id, label=!!prop_sample_labels, prop=!!prop), by="sample_id") %>%
        left_join(samples_colors %>% select(sample_id, clr=!!prop), by="sample_id")
      tib = left_join(tib, samples %>% select(sample_id, exclude), by="sample_id")

      # check if current property (column in sample metadata) only contains numbers. If so, color as numeric
      prop_is_numeric = !prop %in% c("group", "exclude") && suppressWarnings(all(is.finite(as.numeric(tib$prop))))
      if(prop_is_numeric) {
        tib$prop = as.numeric(tib$prop)
      }

      uprop = gtools::mixedsort(unique(tib$prop))
      clr_map = array(paste0(tib$clr[match(uprop, tib$prop)], "CC"), dimnames=list(uprop))

      p = ggplot(tib, aes(x = x, y = y, label = label, colour = prop)) +
        geom_point(aes(shape = ifelse(exclude & pch_as_exclude, "0", "16"))) +
        scale_shape_manual(values = c("0" = 0, "16" = 16), guide = FALSE) +
        labs(
          title = "",
          x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
          y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
        ) +
        theme_bw() #+ theme(plot.title = element_text(hjust = 0.5))

      # recycle colors from upstream calculation
      if(prop_is_numeric) {
        # for gradient colors, sorting by respective value is important (already converted the values to numeric type a few lines above)
        p = p +
          scale_color_gradientn(name = gsub("_", " ", prop), colours = unique(tib$clr[order(tib$prop)])) +
          theme(legend.text = element_text(angle=90, hjust=0.5, vjust=0.5),
                legend.title = element_text(size=10, face = "bold"))
      } else {
        p = p + scale_colour_manual(
          name = gsub("_", " ", prop),
          values = clr_map, # named array, value=color, name=property
          breaks = names(clr_map), # sort the legend
          aesthetics = "colour") +
          guides(colour = guide_legend(title.position = "top")) +
          theme(legend.text = element_text(size = ifelse(length(clr_map) < 6,
                                                         10,
                                                         ifelse(length(clr_map) < 10, 8, 6)) ),
                legend.title = element_text(size=10, face = "bold"))
      }

      p_labeled = p + ggrepel::geom_text_repel(vjust = -.5, show.legend = FALSE, size = 2)

      plotlist[[length(plotlist) + 1]] = p
      plotlist[[length(plotlist) + 1]] = p_labeled
    }

    pcaplots[[length(pcaplots) + 1]] = ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
  }


  ### code snippet for adding PCA on any custom continuous scale. eg, some aspect that reflects sample quality such as; #detect, outliers deviation in retention time, impact on CoV
  # counts = peptides %>% group_by(sample_id) %>% summarise(detect = sum(detect), quant=n()) %>% left_join(samples %>% select(sample_id, shortname), by="sample_id")
  # plotlist = list()
  # for (dims in list(1:2, c(1, 3), 2:3)) { # dims=1:2
  #   tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>% left_join(counts, by="sample_id")
  #
  #   p = ggplot(tib, aes(x = x, y = y, label = shortname, colour = detect)) +
  #     geom_point() +
  #     scale_color_gradient(name = "detected peptides") +
  #     labs(
  #       title = "",
  #       x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
  #       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
  #     ) +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5),
  #           legend.text = element_text(angle=90, hjust=0.5, vjust=0.5))
  #
  #   p_labeled = p + ggrepel::geom_text_repel(show.legend = FALSE, size = 2)
  #
  #   plotlist[[length(plotlist) + 1]] = p
  #   plotlist[[length(plotlist) + 1]] = p_labeled
  # }
  # pcaplots[[length(pcaplots) + 1]] = ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")

  return(pcaplots)
}



#
#
# # for reference, standard pca (does not cope with missing values); summary(stats::prcomp(t(matrix_sample_intensities), center = T, scale. = F))
# plot_sample_pca = function(matrix_sample_intensities, samples, samples_colors, label_samples_by_shortname = TRUE) {
#   PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
#   mat_pca = PPCA@scores # rownames = sample_id
#   pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores)))
#
#   props = grep("^(sample_id$|sample_index$|shortname$|contrast:)", colnames(samples_colors), ignore.case = T, value = T, invert = T)
#
#   # how should we label the samples? if not by shortname, use the index (and fallback to sample_id if index not available)
#   if(length(label_samples_by_shortname) != 1 || !is.logical(label_samples_by_shortname) || !is.finite(label_samples_by_shortname)) {
#     label_samples_by_shortname = nrow(samples) < 20
#   }
#   prop_sample_labels = "shortname"
#   if(!label_samples_by_shortname) {
#     prop_sample_labels = intersect(c("sample_index", "sample_id"), colnames(samples))[1]
#   }
#
#   pcaplots = list()
#   for (prop in props) { #prop=props[4]
#     # don't plot 'exclude' property if there are no excluded samples
#     if(prop == "exclude" && !any(samples$exclude)) {
#       next
#     }
#
#     plotlist = list()
#     for (dims in list(1:2, c(1, 3), 2:3)) { # dims=1:2
#       # extract data from PCA object  &  join with sample metadata to get color-coding and label
#       tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>%
#         left_join(samples %>% select(sample_id, label=!!prop_sample_labels, prop=!!prop), by="sample_id") %>%
#         left_join(samples_colors %>% select(sample_id, clr=!!prop), by="sample_id")
#
#       # check if current property (column in sample metadata) only contains numbers. If so, color as numeric
#       prop_is_numeric = !prop %in% c("group", "exclude") && suppressWarnings(all(is.finite(as.numeric(tib$prop))))
#       if(prop_is_numeric) {
#         tib$prop = as.numeric(tib$prop)
#       }
#
#       uprop = unique(tib$prop)
#       clr_map = array(paste0(tib$clr[match(uprop, tib$prop)], "CC"), dimnames=list(uprop))
#
#       p = ggplot(tib, aes(x = x, y = y, label = label, colour = prop)) +
#         geom_point() +
#         labs(
#           title = "",
#           x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
#           y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
#         ) +
#         theme_bw() #+ theme(plot.title = element_text(hjust = 0.5))
#
#       # recycle colors from upstream calculation
#       if(prop_is_numeric) {
#         # for gradient colors, sorting by respective value is important (already converted the values to numeric type a few lines above)
#         p = p +
#           scale_color_gradientn(name = gsub("_", " ", prop), colours = unique(tib$clr[order(tib$prop)])) +
#           theme(legend.text = element_text(angle=90, hjust=0.5, vjust=0.5),
#                 legend.title = element_text(size=10, face = "bold"))
#       } else {
#         p = p + scale_colour_manual(
#           name = gsub("_", " ", prop),
#           values = clr_map, # named array, value=color, name=property
#           breaks = names(clr_map), # sort the legend
#           aesthetics = "colour") +
#           guides(colour = guide_legend(title.position = "top")) +
#           theme(legend.text = element_text(size = ifelse(length(clr_map) < 6,
#                                                          10,
#                                                          ifelse(length(clr_map) < 10, 8, 6)) ),
#                 legend.title = element_text(size=10, face = "bold"))
#       }
#
#       p_labeled = p + ggrepel::geom_text_repel(vjust = -.5, show.legend = FALSE, size = 2)
#
#       plotlist[[length(plotlist) + 1]] = p
#       plotlist[[length(plotlist) + 1]] = p_labeled
#     }
#
#     pcaplots[[length(pcaplots) + 1]] = ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
#   }
#
#
#   ### code snippet for adding PCA on any custom continuous scale. eg, some aspect that reflects sample quality such as; #detect, outliers deviation in retention time, impact on CoV
#   # counts = peptides %>% group_by(sample_id) %>% summarise(detect = sum(detect), quant=n()) %>% left_join(samples %>% select(sample_id, shortname), by="sample_id")
#   # plotlist = list()
#   # for (dims in list(1:2, c(1, 3), 2:3)) { # dims=1:2
#   #   tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>% left_join(counts, by="sample_id")
#   #
#   #   p = ggplot(tib, aes(x = x, y = y, label = shortname, colour = detect)) +
#   #     geom_point() +
#   #     scale_color_gradient(name = "detected peptides") +
#   #     labs(
#   #       title = "",
#   #       x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
#   #       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
#   #     ) +
#   #     theme_bw() +
#   #     theme(plot.title = element_text(hjust = 0.5),
#   #           legend.text = element_text(angle=90, hjust=0.5, vjust=0.5))
#   #
#   #   p_labeled = p + ggrepel::geom_text_repel(show.legend = FALSE, size = 2)
#   #
#   #   plotlist[[length(plotlist) + 1]] = p
#   #   plotlist[[length(plotlist) + 1]] = p_labeled
#   # }
#   # pcaplots[[length(pcaplots) + 1]] = ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
#
#   return(pcaplots)
# }
