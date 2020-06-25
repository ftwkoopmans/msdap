
#' generate color-coding for all sample metadata
#' @param samples todo
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom gtools mixedsort
#' @importFrom ggplot2 cut_interval
sample_color_coding = function(samples) {
  color_categorical = list(
    # https://observablehq.com/@d3/color-schemes?collection=@d3/d3-scale-chromatic
    # https://personal.sron.nl/~pault/#sec:qualitative  -->> changed order to match recommended for fixed sequence  +  added near-black color at the end (#262626)
    d3_set1 = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf","#999999"), # minus yellow
    d3_dark2 = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666"),

    PaulTol_bright = c(blue="#4477AA",red="#EE6677",green="#228833",yellow="#CCBB44",cyan="#66CCEE",purple="#AA3377",grey="#BBBBBB",darkgrey="#262626"),
    PaulTol_vibrant = c(orange="#EE7733",blue="#0077BB",magneta="#EE3377",cyan="#33BBEE",red="#CC3311",teal="#009988",grey="#BBBBBB", darkgrey="#262626"), # alternative to PaulTol_bright
    # PaulTol_dark = c("#222255","#225555","#225522","#666633","#663333","#555555", darkgrey="#262626"),
    PaulTol_muted_subset = c(rose="#CC6677",indigo="#332288",green="#117733",wine="#882255",teal="#44AA99",olive="#999933",purple="#AA4499", darkgrey="#262626"), # minus cyan and sand

    d3_tabl = c("#4e79a7","#f28e2c","#e15759","#76b7b2","#59a14f","#af7aa1","#9c755f","#bab0ab"), # minus light yellow and light pink
    d3_cat10 = rev(c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f")) # minus last 2 ,"#bcbd22"
  )
  # # debug; show colors;
  # for(i in seq_along(color_categorical)) {
  #   barplot(rep(1,length(color_categorical[[i]])), col=color_categorical[[i]], border = NA, yaxt='n', main=names(color_categorical)[i])
  #   plot(runif(length(color_categorical[[i]])*3), runif(length(color_categorical[[i]])*3), col=rep(color_categorical[[i]], 3), pch=16, xaxt="n", yaxt="n", xlab="", ylab="", main=names(color_categorical)[i])
  # }

  # d3_set2 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854", "#e5c494","#b3b3b3"), # minus yellow -->> too flat
  #"Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999"), # RColorBrewer::brewer.pal(9, "Set1")   minus yellow
  #"Set2" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"), # RColorBrewer::brewer.pal(8, "Set2")   as-is
  #"Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), # RColorBrewer::brewer.pal(8, "Dark2")   as-is
  #"Accent" = c("#666666", "#BF5B17", "#F0027F", "#386CB0", "#FFFF99", "#FDC086", "#BEAED4", "#7FC97F"), # rev(RColorBrewer::brewer.pal(8, "Accent"))


  #
  # color_pals = tibble(name = c("Set1", "Set2", "Dark2", "Accent"), reverse = c(F, F, F, T))
  color_pals_continuous = tibble(name = c("GnBu", "YlOrRd", "YlGn")) # c("BuPu", "OrRd", "YlGn")

  # color_pals$maxcolors = RColorBrewer::brewer.pal.info$maxcolors[match(color_pals$name, rownames(RColorBrewer::brewer.pal.info))]
  # if(any(is.na(color_pals$maxcolors))) {
  #   append_log(paste("unknown color palette;", paste(color_pals$name[is.na(color_pals$maxcolors)], collapse = ", ")), type = "error")
  # }

  ### color palette per property
  sample_colors = samples %>% select(sample_id, shortname)
  sample_properties = grep("^(sample_id$|sample_index$|shortname$|condition$|(all_){0,1}_proteins$|^key_|contrast:)", colnames(samples), value = T, ignore.case = T, invert = T)

  count = 0
  count_continuous = 0
  for (i in seq_along(sample_properties)) {
    prop = sample_properties[i]
    # values to color-code
    val = samples %>% select(!!prop) %>% pull()
    uval = unique(val)

    # if not color-coding group/exclude, we skip if there is nothing to color-code
    if(length(uval) == 1 && !prop %in% c("group", "exclude")) {
      next
    }

    prop_is_numeric = !prop %in% c("group", "exclude") && suppressWarnings(all(is.finite(as.numeric(uval))))
    if(prop_is_numeric) {
      val = as.numeric(val)
      ## create continuous color gradient between lowest and highest value
      # cycle through color palettes
      palette_index = 1 + (count_continuous %% nrow(color_pals_continuous))
      count_continuous = count_continuous + 1
      # standard color scale, but take out 2 colors closest to white
      palette_colors = RColorBrewer::brewer.pal(9, color_pals_continuous$name[palette_index])[-1:-2]

      if(length(uval) == 1) {
        # if there is just one unique value, use the last (darkest) value in this palette
        sample_colors[,prop] = tail(palette_colors, 1)
      } else {
        # interpolate to more refined
        clr = colorRampPalette(palette_colors)(25)
        # optionally, limit outliers/extremes
        q = quantile(val, probs = c(.01, .99))
        val[val < q[1]] = q[1]
        val[val > q[2]] = q[2]
        # map values to colors
        val_cut = as.numeric(ggplot2::cut_interval(val, length(clr)))
        sample_colors[,prop] = clr[val_cut]
      }
    } else {
      if(prop != "group") {
        uval = gtools::mixedsort(uval)
      }
      # cycle through color palettes
      palette_index = 1 + (count %% length(color_categorical)) #(count %% nrow(color_pals))
      count = count + 1
      # create a color palette, optionally interpolating if there are more unique values than palette size
      uval_clr = color_categorical[[palette_index]]
      if(length(uval) > length(uval_clr)) {
        uval_clr = colorRampPalette(uval_clr)(length(uval))
      }

      # store respective colors
      sample_colors[,prop] = uval_clr[match(val, uval)]
    }

  }
  return(sample_colors)
}



#' convert color matrix to long format
#' @param samples todo
#' @param samples_colors todo
sample_color_coding_as_long_tibble = function(samples, samples_colors) {
  # iterate each property used for color coding, then combine the property, value and respective color in long-format tibble
  tib_color_coding = NULL
  for (prop in setdiff(colnames(samples_colors), c("sample_id", "sample_index", "shortname"))) { #prop="group"
    tib_color_coding = bind_rows(tib_color_coding, samples %>%
                                   select(shortname, val=!!prop) %>%
                                   mutate_all(as.character) %>%
                                   add_column(clr = samples_colors %>% pull(!!prop)) %>%
                                   add_column(prop = prop, .after = "shortname"))
  }

  # exclude color codings that have just 1 unique value
  tib_counts = tib_color_coding %>% group_by(prop) %>% summarise(n=n_distinct(clr))
  tib_color_coding = tib_color_coding %>% filter(prop %in% (tib_counts %>% filter(n>1) %>% pull(prop)))

  # convert the properties to a factor, ordered exactly as in the input sample metadata table (eg; respective column names)
  tib_color_coding$prop = factor(tib_color_coding$prop, levels = intersect(colnames(samples), tib_color_coding$prop))

  return(tib_color_coding)
}



#' placeholder title
#' @param samples todo
#' @param samples_colors_long todo
#'
#' @importFrom ggpubr theme_pubr
ggplot_sample_detect_counts_barplots = function(samples, samples_colors_long) {
  tib = samples

  # if there are 'all peptide' counts (detect + MBR/quant-only), include these in plot tibble
  if("all_peptides" %in% colnames(samples)) {
    tib$quantified_peptides = tib$all_peptides - tib$detected_peptides
    tib$quantified_proteins = tib$all_proteins - tib$detected_proteins
  }

  # from wide to long format
  tib = tib %>%
    select(!!c("shortname", grep("(detected|quantified)_(peptides|proteins)", colnames(tib), value = T))) %>%
    pivot_longer(cols = -shortname, names_to = "type", values_to = "count")

  # flip levels/sorting because we use coord_flip() downstream in ggplot
  tib$shortname = factor(tib$shortname, levels = rev(samples$shortname))
  tib$pep_or_prot = ifelse(grepl("peptide", tib$type), "peptides", "proteins")
  tib$detect_or_quant = ifelse(grepl("detected", tib$type), "detect", "quant")

  # rename labels for plot clarity
  tib$type["detected_peptides" == tib$type] = "peptides: identified & quantified"
  tib$type["quantified_peptides" == tib$type] = "peptides: only quantified"
  tib$type["detected_proteins" == tib$type] = "proteins: identified & quantified"
  tib$type["quantified_proteins" == tib$type] = "proteins: only quantified"

  # color-coding
  # clr = c("detected_peptides"="#0570b0", "quantified_peptides"="#74a9cf", "detected_proteins"="#6a51a3", "quantified_proteins"="#9e9ac8")
  clr = c("peptides: identified & quantified"="#0570b0", "peptides: only quantified"="#74a9cf",
          "proteins: identified & quantified"="#6a51a3", "proteins: only quantified"="#9e9ac8")

  # coordinates for our customized dual-barplot
  tib = tib %>% group_by(pep_or_prot, shortname) %>% mutate(total=sum(count)) %>% arrange(shortname, type)
  tib_dual = tib %>% group_by(pep_or_prot) %>% mutate(count_scaled = count / max(total) * ifelse(pep_or_prot=="proteins", -1, 1),
                                                      count_scaled_max = total / max(total) * ifelse(pep_or_prot=="proteins", -1, 1))
  tib_dual$outlier_lowside = abs(tib_dual$count_scaled_max) < 0.3

  # custom x-axis labels, since this dimension not ranges from -1:1
  ticks_peptides = pretty(1:max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total)), n = 5, min.n = 4)
  ticks_peptides_scaled = ticks_peptides / max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total))
  ticks_proteins = pretty(1:max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total)), n = 5, min.n = 4)
  ticks_proteins_scaled = ticks_proteins / max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total))
  # combine
  ticks_coord = c(-1 * rev(ticks_proteins_scaled), ticks_peptides_scaled[-1])
  ticks_label = c(rev(ticks_proteins), ticks_peptides[-1])


  # plot code is somewhat convoluted by careful alignment of text labels. simplest QC plot of data as-is: ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) + geom_bar(stat = "identity", position = position_stack(reverse = T)) + coord_flip()
  p = ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T)) +
    geom_text(aes(label = count,
                  y = ifelse(detect_or_quant=="detect", ifelse(pep_or_prot=="proteins", -0.025, 0.025), count_scaled_max),
                  hjust = ifelse((type %in% c("peptides: identified & quantified", "proteins: only quantified") & !c(type=="proteins: only quantified" & outlier_lowside)) |
                                   (type %in% c("peptides: only quantified") & outlier_lowside), -0.1, 1.1)),
              colour = ifelse(tib_dual$outlier_lowside & tib_dual$detect_or_quant=="quant", "darkgrey", "white"),
              size = 3.5, # 4 is default
              check_overlap = F) +
    scale_y_continuous(breaks=ticks_coord, labels=ticks_label, ) +
    scale_fill_manual(values = clr, guide = guide_legend(nrow = 2)) +
    coord_flip() +
    labs(x = "", y = "protein / peptide counts (axis is mirrored, with separate scales for each") +
    ggpubr::theme_pubr(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(0.75,"line"), legend.text = element_text(size=9),
          panel.grid = element_blank(), axis.line.y.right = element_blank(),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

  return(p)
  ###### some reference code for separate plots
  # tib_pep = tib %>% filter(pep_or_prot == "peptides") %>% droplevels()
  # tib_prot = tib %>% filter(pep_or_prot == "proteins") %>% droplevels()
  #
  # plot_text_lim = 0.1 * max(tib_pep$count)
  # p_pep = ggplot(tib_pep, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
  #
  # plot_text_lim = 0.1 * max(tib_prot$count)
  # p_prot = ggplot(tib_prot, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
}



#' compact scatterplot with color-coding captures all sample metadata in one figure
#' @param tib_plot todo
ggplot_sample_detect_vs_metadata_scatterplot = function(tib_plot) {
  # by appending the property to respective values we ensure that duplicated values across the metadata table (eg; 2 properties that both take TRUE/FALSE values) are color-coded properly
  tib = tib_plot %>% mutate(val = paste(prop,val))
  # unique set of colors used in this plot
  tib_col = tib %>% distinct(val, clr)

  # reverse order of properties
  levels(tib$prop) <- rev(levels(tib$prop))

  # pseudocode illustrating how we revert y-axis indice order; y = 1 + (1:4 - 4) + jitter
  tib$prop_y = 1 + abs(as.numeric(tib$prop) - length(levels(tib$prop))) + tib$sample_jitter

  return(ggplot(tib, aes(x=detected_peptides, y=prop, colour = val, fill = val)) +
           geom_point(alpha = 0) + # first plot invisible points to fix y axis
           geom_point(aes(y=prop_y), shape = 21) + # then plot with custom y location for our fixed jitter
           scale_colour_manual(values = array(paste0(tib_col$clr, "BB"), dimnames = list(tib_col$val)), aesthetics = c("colour")) +
           scale_colour_manual(values = array(paste0(tib_col$clr, "66"), dimnames = list(tib_col$val)), aesthetics = c("fill")) +
           theme_bw() +
           labs(x="number of detected peptides per sample", y="", title="detection rates per sample color-coded by sample metadata") +
           theme(legend.position = "none", plot.title = element_text(size = 10), axis.title.x.bottom = element_text(size = 9)))

  ### example code, externalized upstream
  # # problem: default ggplot jitter is not recycled (eg; same datapoint in multiple groups has distinct jitter, making them hard to visually align/compare)
  # # solution: precalculate jitter
  # # 1) some jitter for each sample
  # tib$shortname = as.factor(tib$shortname)
  # jit = runif(seq_along(levels(tib$shortname)), -0.2, 0.2)
  # # 2) map to plot tibble; y location = index of property + jitter of respective sample
  # # pseudocode illustrating how we revert y-axis indice order; y = 1 + (1:4 - 4) + jitter
  # tib$prop_y = 1 + abs(as.numeric(tib$prop) - length(levels(tib$prop))) + jit[as.numeric(tib$shortname)]
  # # tib$prop_y = as.numeric(tib$prop) + jit[as.numeric(tib$shortname)] # default, no re-ordering of properties

}



#' separate plot for each metadata property
#' @param tib_plot todo
ggplot_sample_detect_vs_metadata_scatterplot_by_prop = function(tib_plot) {
  result = list()

  for(prp in levels(tib_plot$prop)) { # prp="group"
    # subset plot data for current metadata property
    tib = tib_plot %>% filter(prop == prp)
    # median detect count by unique value
    tib_stat = tib %>%
      group_by(val) %>%
      summarise(median_detect_count=median(detected_peptides)) %>%
      arrange(desc(median_detect_count))
    tib_stat$index = 1:nrow(tib_stat)

    # merge counts with plot tibble and convert 'val' to factor to enforce sorting
    tib = left_join(tib, tib_stat, by="val")
    tib$val = factor(tib$val, levels = tib_stat$val)
    # unique set of colors used in this plot
    tib_col = tib %>% distinct(val, clr)

    # see ggplot_sample_detect_vs_metadata_scatterplot()
    tib$prop_y = as.numeric(tib$val) + tib$sample_jitter

    result[[prp]] = list(n = nrow(tib_stat), # n = number of rows/elements. useful for deciding on plot size downstream
                         plot = ggplot(tib, aes(x=detected_peptides, y=val, colour = val, fill = val)) +
                           geom_point(alpha = 0) + # first plot invisible points to fix y axis
                           geom_point(aes(y=prop_y), shape = 21) + # then plot with custom y location for our fixed jitter
                           geom_segment(data=tib_stat, aes(x=median_detect_count, xend=median_detect_count, y=index+0.3, yend=index-0.3, colour=val)) + # x = 2, y = 15, xend = 3, yend = 15
                           scale_colour_manual(values = array(paste0(tib_col$clr, "BB"), dimnames = list(tib_col$val)), aesthetics = c("colour")) +
                           scale_colour_manual(values = array(paste0(tib_col$clr, "66"), dimnames = list(tib_col$val)), aesthetics = c("fill")) +
                           # scale_colour_manual(values = array(tib_col$clr, dimnames = list(tib_col$val)), aesthetics = c("colour", "fill")) +
                           theme_bw() +
                           labs(x="number of detected peptides per sample", y="", title=paste("sample metadata used for color-coding:", prp)) +
                           theme(legend.position = "none", plot.title = element_text(size = 10), axis.title.x.bottom = element_text(size = 9)) )
  }

  return(result)
}
