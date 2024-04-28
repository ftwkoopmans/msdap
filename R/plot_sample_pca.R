
#' placeholder title
#' for reference, standard pca (does not cope with missing values); summary(stats::prcomp(t(matrix_sample_intensities), center = T, scale. = F))
#' @param matrix_sample_intensities peptide or protein abundance matrix
#' @param samples sample metadata table
#' @param samples_colors sample colors in wide format tibble
#' @param sample_label_property labels used to represent the samples. Set to "auto" to automatically select based on number of samples (default). possible values; auto, shortname, index, index_asis
#' @param pch_as_exclude should exclude samples be represented with a distinct symbol shape? possible values; FALSE, TRUE (default)
#' @param infer_continuous_scale should the legends infer continuous scale from the data? FALSE = all categorical, TRUE = automatic (default)
#' @importFrom pcaMethods pca
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom gtools mixedsort
#' @importFrom scales rescale
plot_sample_pca = function(matrix_sample_intensities, samples, samples_colors, sample_label_property = "auto", pch_as_exclude = TRUE, infer_continuous_scale = TRUE) {
  if(!sample_label_property %in% c("auto", "shortname", "index", "index_asis")) {
    append_log(paste("invalid value for parameter 'sample_label_property':", sample_label_property), type = "error")
  }
  if(!all(c("sample_index", "sample_id", "shortname") %in% colnames(samples))) {
    append_log(paste(paste(c("sample_index", "sample_id", "shortname"), collapse = ", "), "are required column in the samples table"), type = "error")
  }

  PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
  mat_pca = PPCA@scores # rownames = sample_id
  pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores)))

  props = grep("^(sample_id$|sample_index$|shortname$|contrast:)", colnames(samples_colors), ignore.case = T, value = T, invert = T)

  # if more than 50 samples, don't ggrepel the sample labels and instead simply plot each label at the PCA x/y location (color-coded by metadata as per usual). No point/shape is drawn in this case
  textonly_sample_labels = nrow(samples) > 50 # default
  prop_sample_labels = "shortname" # default
  if(sample_label_property == "auto") {
    # only shortname if less than 25 samples AND 80% of those are actually short strings (<=30 characters)
    prop_sample_labels = ifelse(nrow(samples) < 25 && sum(nchar(samples$shortname) > 30) < nrow(samples) * 0.2, "shortname", "sample_index") # alternatively; nrow(samples) * stats::median(nchar(samples$shortname)) < 25*10
  }
  if(sample_label_property == "index") {
    prop_sample_labels = "sample_index"
  }
  if(sample_label_property == "index_asis") {
    prop_sample_labels = "sample_index"
    textonly_sample_labels = TRUE
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

      plot_as_numeric = FALSE
      if(infer_continuous_scale) {
        # check if current property (column in sample metadata) only contains numbers
        prop_is_numeric = suppressWarnings(all(is.finite(as.numeric(na.omit(tib$prop)))))
        # make a numeric plot only if there are more than 3 unique values. Also covers the case of booleans that should be plotted as categorical variables
        plot_as_numeric = !prop %in% c("group", "exclude") && prop_is_numeric && n_distinct(as.numeric(tib$prop), na.rm = T) > 2
      }

      # mutate the actual data before making ggplot object
      if(plot_as_numeric) {
        tib$prop = as.numeric(tib$prop)
      } else {
        # for categorical variables, convert NA values to a character so they show up in the legend
        tib$prop[is.na(tib$prop)] = "<NA>"
        tib$prop = as.character(tib$prop)
      }

      # base plot
      p = ggplot(tib, aes(x = x, y = y, label = label, colour = prop, fill = prop)) +
        guides(alpha = "none", fill = "none") +
        labs(
          title = "",
          x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
          y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
        ) +
        # coord_cartesian(clip = "off") +
        theme_bw()


      if(plot_as_numeric) {
        tib = tib %>%
          mutate(
            # add alpha to colors
            clr_fill = paste0(clr, "66"),
            clr = paste0(clr, "BB"),
            # scale color values between 0~1
            prop = scales::rescale(prop)) %>%
          # sort values from low to high
          # doesn't matter that the order of the values submitted to ggplot2::scale_color_X is changed from when we created ggplot object above
          arrange(prop)

        clr_na = "#BBBBBBBB"
        if(any(is.na(tib$prop))) {
          # get NA color from data
          clr_na = tib %>% filter(is.na(prop)) %>% pull(clr) %>% head(1)
          # set color code to NA or they'll show up in ggplot's color scale / legend
          tib$clr[is.na(tib$prop)] = NA
          tib$clr_fill[is.na(tib$prop)] = NA
        }


        # for gradient colors, sorting by respective value is important (already converted the values to numeric type a few lines above)
        p = p +
          scale_color_gradientn(name = gsub("[ _]+", " ", prop),
                                values = tib$prop, # set the values parameter to use the exact value-to-color mapping we computed upstream
                                colours = tib$clr,
                                aesthetics = "colour",
                                na.value = clr_na) +
          scale_color_gradientn(values = tib$prop,
                                colours = tib$clr_fill,
                                aesthetics = "fill",
                                na.value = clr_na) +
          theme(legend.text = element_text(angle=90, hjust=0.5, vjust=0.5),
                legend.title = element_text(size=10, face = "bold"))
      } else {
        uprop = gtools::mixedsort(unique(tib$prop))
        clr_map = array(tib$clr[match(uprop, tib$prop)], dimnames=list(uprop))

        p = p +
          scale_colour_manual(
            name = gsub("[ _]+", " ", prop),
            values = paste0(clr_map, "BB"), # named array, value=color, name=property
            breaks = names(clr_map), # sort the legend
            aesthetics = "colour") +
          scale_colour_manual(
            values = paste0(clr_map, "66"), # named array, value=color, name=property
            breaks = names(clr_map), # sort the legend
            aesthetics = "fill") +
          guides(colour = guide_legend(title.position = "top")) +
          theme(legend.text = element_text(size = ifelse(length(clr_map) < 6,
                                                         10,
                                                         ifelse(length(clr_map) < 10, 8, 6)) ),
                legend.title = element_text(size=10, face = "bold"))
      }

      if(textonly_sample_labels) {
        p_labeled = p + geom_text(show.legend = FALSE, size = 2)
      } else {
        p_labeled = p +
          geom_point(aes(shape = I(ifelse(exclude & pch_as_exclude, 0, 21)))) +
          ggrepel::geom_text_repel(show.legend = FALSE, size = 2, segment.alpha = .3, min.segment.length = 0, na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, box.padding = 0.2, seed = 123)
      }

      plotlist[[length(plotlist) + 1]] = p + geom_point(aes(shape = I(ifelse(exclude & pch_as_exclude, 0, 21))))
      plotlist[[length(plotlist) + 1]] = p_labeled
    }

    # finally, collapse individual PCA plots
    # importantly: if samples are labeled by index, add a warning/note
    # eg; if user labels samples 1~6 as shortname and our dataset$samples table has a different order (causing sample_index to not align), users may mistake sample identities
    # users labeling samples with shortnames like  WT:1,2,4 KO:5,6,7  happens but all too often, so this note must be visible on same page as the plot
    plotlist_merged = ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
    if(prop_sample_labels == "sample_index") {
      plotlist_merged = ggpubr::annotate_figure(plotlist_merged, top = ggpubr::text_grob('samples are labeled by "sample_index" (described in samples.xlsx included with MS-DAP output)', size = 10))
    }
    pcaplots[[prop]] = plotlist_merged
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





#' helper function for use-case beyond making a report using standard analysis_quickstart() workflow; PCA on the subset of samples-of-interest from some contrast
#'
#' @param dataset your dataset after application of dea. eg; results from analysis_quickstart()
#' @param contr the exact name of the contrast. Example; "contrast: ctrl1,ctrl2 vs phenotype1,phenotype2,phenotype3"  This should be one of the columns in dataset$samples table (which are generated by the setup_contrasts() function)
#' @param sample_label_property see plot_sample_pca() function for params
#' @importFrom tidyr pivot_wider
#'
#' @export
plot_sample_pca__sample_in_contrast = function(dataset, contr, sample_label_property = "auto") {
  if(length(contr) != 1 || !(contr %in% colnames(dataset$peptides) || contr %in% colnames(dataset$samples)) ) {
    append_log("The contrast has to be a column name in the dataset's samples table, as created by the setup_contrasts function. Typical use-case: setup_contrasts(), analysis_quickstart(), then call this function.", type="error")
  }

  if(contr %in% colnames(dataset$peptides)) {
    intensity_col_contr = contr
  } else {
    intensity_col_contr = paste0("intensity_", contr)
    if(! intensity_col_contr %in% colnames(dataset$peptides)) {
      append_log(sprintf("The column '%s' cannot be found in the dataset's peptides table, by_contrast filtering has to be enabled for this function to work (it relies on the subsetting of samples + filtering + normalization procedures that create a data matrix specific to some statistical contrast). Typical use-case: setup_contrasts(), analysis_quickstart(), then call this function.", intensity_col_contr), type="error")
    }
  }

  # sample color-coding, exactly like report_as_rmarkdown.R  (should compute colors based on entire dataset so color-coding is the same as the report that features the full dataset)
  samples_colors_long = sample_color_coding__long_format(dataset$samples)
  samples_colors = samples_colors_long %>% dplyr::select(sample_id, shortname, prop, clr) %>% tidyr::pivot_wider(id_cols = c(sample_id, shortname), names_from=prop, values_from=clr)
  # case data to wide-format, then plot PCA
  tibw = dataset$peptides %>%
    dplyr::select(key_peptide, sample_id, intensity=!!intensity_col_contr) %>% dplyr::filter(!is.na(intensity)) %>%
    tidyr::pivot_wider(id_cols = key_peptide, names_from = sample_id, values_from = intensity)
  # return the results of the respective plot function
  return( suppressWarnings(plot_sample_pca(as_matrix_except_first_column(tibw), samples = dataset$samples, samples_colors = samples_colors, sample_label_property = sample_label_property)) )
}



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
