
#' placeholder title
#' @param tib_input todo
#' @param samples todo
#' @param isdia todo
plot_abundance_distributions = function(tib_input, samples, isdia) {
  param_density_bandwidth = "sj"
  param_density_adjust = 1.0 # alternatively, 0.9

  # for DIA datasets, only use abundances for peptides detected at some q-value threshold (instead of 'mbr' at any confidence threshold)
  tib = left_join(tib_input %>% filter(detect | !isdia),
                  samples %>% select(sample_id, shortname, group, exclude),
                  by = "sample_id") %>%
    mutate(intensity = intensity / log2(10)) # change log base from log2 to log10
  overall_xlim = quantile(tib$intensity, probs = c(.001, .999), na.rm = T)

  tib$exclude = factor(tib$exclude, levels = c(FALSE, TRUE))

  ## all samples in a single plot
  p_all_intensity_distributions = ggplot(tib, aes(x=intensity, colour=shortname, labels=shortname, linetype = exclude)) +
    geom_line(stat = "density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm = T) +
    scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
    coord_cartesian(xlim = overall_xlim) +
    labs(x="peptide intensity (log10)", y="density") +
    facet_grid(~"all samples") +
    theme_bw() +
    theme(legend.position = "none", legend.title = element_blank())

  # overall distributions, per group so plots are readable for large datasets
  l_group_intensity_distributions = list()
  for (g in unique(samples$group)) {
    n_samples = sum(samples$group == g)
    p = ggplot(tib %>% filter(group == g), aes(x = intensity, colour = shortname, labels = shortname, linetype = exclude)) +
      stat_density(geom="line", position="identity", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm = T) +
      scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
      coord_cartesian(xlim = overall_xlim) +
      labs(x="peptide intensity (log10)", y="density", colour="sample", linetype = "exclude sample") +
      guides(linetype = "none", #guide_legend(direction = "horizontal", title.position = "top"),
             colour = guide_legend(title = NULL)) +
      facet_grid(~group) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = ifelse(n_samples < 4, 10, ifelse(n_samples < 6, 8, 6)) ),
            legend.title = element_text(size=10) )

    ## split legend into separate plot if more than N entries
    if(n_samples > 12) {
      p_split = ggplot_split_legend(p)
      l_group_intensity_distributions[[g]] = p_split$plot
      l_group_intensity_distributions[[paste0(g, "_legend")]] = p_split$legend # + guides(colour = guide_legend(ncol = 2))
    } else {
      l_group_intensity_distributions[[g]] = p
    }

  }

  return(list(intensity_distributions_all = p_all_intensity_distributions, intensity_distributions_bygroup = l_group_intensity_distributions))
}




#' Plot within-group foldchange distributions to assess outlier samples
#'
#' @examples
#' \dontrun{
#'   dataset = filter_peptides_by_group(
#'     dataset,
#'     colname="intensity_qc_basic",
#'     disregard_exclude_samples = FALSE,
#'     nquant=2,
#'     ndetect=1+is_dia_dataset(dataset),
#'     norm_algorithm = c("vsn", "modebetween_protein"),
#'     rollup_algorithm = "maxlfq"
#'   )
#'   l = plot_foldchange_distribution_among_replicates(mypeptides, dataset$samples)
#'   for(p in l$plotlist) print(p)
#'   print(l$tib_scores)
#' }
#' @param peptides peptide tibble, e.g. dataset$peptides, which must hold raw intensity data in "intensity" column and filtered and normalized data in the column "intensity_qc_basic"
#' @param samples sample tibble, e.g. dataset$samples
plot_foldchange_distribution_among_replicates = function(peptides, samples) {
  param_density_bandwidth = "sj"
  param_density_adjust = 1.0 # alternatively, 0.9

  plotlist = list()
  tib_scores = NULL
  ugrp = unique(samples$group)
  for (g in ugrp) { # g = ugrp[1]
    g_samples = samples %>% filter(group == g)
    skey = g_samples$key_sample
    # cannot make within-group foldchange plots for groups with less than 2 non-exclude samples
    if(sum(g_samples$exclude == FALSE) < 2) {
      append_log(sprintf("skipping within-group foldchange plots for sample group '%s', require at least 2 replicates that are not flagged as 'exclude'", g), type = "info")
      next
    }

    # if a peptide*sample data point has a filtered+normalized intensity, it has at least N observations within this group
    # this filtering step ensures we only data points that are in both the input and the output (eg; same filtering applied to both, just normalization that differs)
    DT = data.table::setDT(peptides %>% select(key_peptide, key_sample, intensity, intensity_qc_basic) %>% filter(key_sample %in% skey & is.finite(intensity_qc_basic)))

    # cannot make within-group foldchange plots for groups that contain only exclude samples, or where all peptide intensities were set to "NA" following filtering
    if(n_distinct(DT$key_peptide) < 50) { # require at least 50 unique peptides
      append_log(sprintf("skipping within-group foldchange plots for sample group '%s', not enough data (from non-exclude samples) to create plots", g), type = "info")
      next
    }

    ## subtract the mean value per peptide from each peptide*sample data point
    ## update for MS-DAP 1.7; now using tmean instead of mean
    # more robust estimates by using trimmed mean (when there are many samples)
    # alternatively, use a robust estimator like MASS::hubers (when sufficient replicates are available)
    # median is not a good choice, we observed cases of strong bias in various small test datasets (e.g. 2 or 3 reps)
    DT_diff_to_mean = DT[ , .(key_sample=key_sample,
                              input = intensity - mean(intensity, na.rm=T, trim=0.2),
                              normalized = intensity_qc_basic - mean(intensity_qc_basic, na.rm=T, trim=0.2)),
                          by=key_peptide]

    # convert data we efficiently gathered with data.table into a long-format tibble
    tib_plot = as_tibble(data.table::melt.data.table(DT_diff_to_mean, id.vars = c("key_peptide", "key_sample"), variable.name = "normalization_status", value.name = "logratio")) %>%
      # join sample metadata (shortname and exclude)
      left_join(samples %>% select(key_sample, shortname, exclude), by = "key_sample") %>%
      # convert 'exclude' into a factor with fixed levels
      mutate(exclude = factor(exclude, levels = c(FALSE, TRUE)),
             group = g) # we already have group variable so we can simply assign it here, don't need to collect in the above left-join (although that would be fine as well, performance difference is moot in this use-case)

    # plot limits
    g_fc_xlim = c(-1, 1) * max(abs(quantile(DT_diff_to_mean$input, probs = c(.005, .995), na.rm = T)))

    ### v2; updated/extended plots (monochrome + highlight topN outliers + density)

    # color-code for highlighted samples
    # number of colors in this array dictates the topN samples to be highlighted
    # color code array, from relatively 'best' to 'worst' sample
    clr_ref = tibble(
      clr = c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#FA9FB5", "#D86FEF", "#6A3D9A", "#FDBF6F", "#FF7F00", "#E31A1C"),
      # if only N colors are in the palette, which should we use?
      prio = c(5, 3, 6, 2, 9, 4, 8, 10, 1, 7)) %>%
      # the ordering of colors; from 'worst' to 'best' sample
      mutate(order = rev(dplyr::row_number()))
    #debug; scales::show_col(clr_ref$clr); scales::show_col(clr_ref %>% filter(prio <= 4) %>% pull(clr)); scales::show_col(clr_ref %>% filter(prio <= 6) %>% pull(clr)); scales::show_col(clr_ref %>% filter(prio <= 8) %>% pull(clr))
    clr_default = "#888888"

    # base color-codes on an outlier metric: standard deviation of normalized intensities
    tib_sample_score = tib_plot %>%
      filter(normalization_status == "normalized") %>%
      group_by(shortname) %>%
      summarise(sd = sd(logratio[logratio > quantile(logratio, probs=0.005) & logratio < quantile(logratio, probs=0.995)]), .groups = "drop") %>%
      # summarise(sd = sd(logratio), .groups = "drop") %>%
      # sort such that 'worst' samples are at the bottom
      arrange(sd) %>%
      # flag bottom N samples
      mutate(has_color = dplyr::row_number()  >  n() - nrow(clr_ref),
             color = clr_default)

    tib_scores = bind_rows(tib_scores, tib_sample_score %>% mutate(group = g) %>% arrange(desc(sd)))

    # assign color to samples that need to be color-coded
    # note that this code works in tandem with above 'has_color' definition; the number of flagged samples is always <= number of available colors in clr_ref
    ncolor = sum(tib_sample_score$has_color)
    tib_sample_score$color[tib_sample_score$has_color] = clr_ref %>% filter(prio <= ncolor) %>% pull(clr)

    # enforce sorting of samples in the plot by setting factor levels
    tib_plot = tib_plot %>%
      left_join(tib_sample_score, by = "shortname") %>%
      mutate(shortname = factor(shortname, levels = tib_sample_score$shortname))


    # monochrome plot without legend
    p_mono = ggplot(tib_plot, aes(x = logratio, labels = shortname, linetype = exclude)) +
      geom_vline(xintercept = 0, size = 1, colour = "darkgrey") +
      geom_line(stat = "density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm = T, alpha = 0.33, size = .75) +
      scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
      facet_wrap(.~group + normalization_status) +
      xlim(g_fc_xlim) +
      labs(x = "log2 fold-change compared to group mean", colour = "sample", linetype = "exclude sample") +
      theme_bw() +
      theme(legend.position = "none")


    # color-code the topN 'worst' samples and show these in the legend
    p_highlight = ggplot(tib_plot, aes(x = logratio, labels = shortname, colour = shortname, linetype = exclude)) +
      geom_vline(xintercept = 0, size = 1, colour = "darkgrey") +
      # first, plot the non-outlier samples (color black, strong transparency)
      geom_line(data = tib_plot %>% filter(has_color == F), stat = "density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm = T, alpha = 0.33, size = 0.75) +
      # next, draw the highlighted samples (aka. worst sd()) on top. Note that these were sorted upstream to ensure the 'worst' sample is drawn on top. (color-coded lines, no transparency)
      geom_line(data = tib_plot %>% filter(has_color == T), stat = "density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm = T, alpha = 1, size = 1) +
      scale_colour_manual(values = array(tib_sample_score$color, dimnames = list(tib_sample_score$shortname)),
                          # the sample labels to be shown in the legend = highlighted samples from summary statistic tibble, in reverse order (so the 'worst' sample comes first)
                          breaks = tib_sample_score %>% filter(has_color==T) %>% pull(shortname) %>% rev() ) +
      scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
      facet_wrap(.~group + normalization_status) +
      xlim(g_fc_xlim) +
      labs(x = "log2 fold-change compared to group mean", colour = "sample", linetype = "exclude sample") +
      guides(linetype = "none", colour = guide_legend(title = element_blank(), ncol = 2, byrow = F)) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7))

    # split legend from the plot, such that both plots have the same dimensions
    p_split = ggplot_split_legend(p_highlight)

    plotlist[[g]] = list(
      mono = p_mono,
      highlight = p_split$plot,
      legend = p_split$legend
    )


    ### v1
    # p = ggplot(tib_plot, aes(x = logratio, labels = shortname, colour = shortname, linetype = exclude)) +
    #   geom_vline(xintercept = 0, size = 1, colour = "darkgrey") + # , linetype = 2
    #   geom_line(stat = "density", size = 1, alpha=0.8, show.legend = T, na.rm=T) +
    #   scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
    #   facet_wrap(.~group + normalization_status) +
    #   xlim(g_fc_xlim) +
    #   labs(x = "log2 fold-change compared to group mean", colour = "sample", linetype = "exclude sample") +
    #   guides(linetype = "none", colour = guide_legend(title = element_blank())) +
    #   theme_bw() +
    #   theme(legend.position = "bottom",
    #         legend.text = element_text(size = ifelse(length(skey) < 4, 10, ifelse(length(skey) < 6, 8, 6)) ),
    #         legend.title = element_text(size=10))
    #
    # ## split legend into separate plot if more than N entries
    # if(dplyr::n_distinct(tib_plot$shortname) > 12) {
    #   p_split = ggplot_split_legend(p)
    #   plotlist[[g]] = p_split$plot
    #   plotlist[[paste0(g, "_legend")]] = p_split$legend # + guides(colour = guide_legend(ncol = 2))
    # } else {
    #   plotlist[[g]] = p
    # }
  }

  return(list(plotlist = plotlist, tib_scores = tib_scores))
}



#' placeholder title
#' @param tib_input todo
#' @param samples todo
#' @param samples_colors todo
#'
#' @importFrom tibble enframe
#' @importFrom ggrepel geom_text_repel
#' @importFrom matrixStats rowSums2 colMedians
ggplot_coefficient_of_variation__leave_one_out = function(tib_input, samples, samples_colors) {
  param_density_bandwidth = "sj"
  param_density_adjust = 1.0 # alternatively, 0.9

  start_time = Sys.time()
  # XX<<-tib_input; YY<<-samples; ZZ<<-samples_colors
  # tib_input=XX; samples=YY; samples_colors=ZZ

  ### leave-one-out computation, per group
  tib_loo_cov = tibble()
  plotlist = list()
  for(grp in unique(samples$group)) { #grp=unique(samples$group)[1]
    # samples for current group
    sid = samples %>% filter(group == grp) %>% pull(sample_id)
    # we use the basic filter+normalization data as input, as this has already been pre-selected for bare minimum of detect/quant (less so than filters below) AND normalized, making CoV estimations more reflective of actual data
    # tibw_grp_intensities = tib_input %>% filter(sample_id %in% sid & !is.na(intensity_qc_cov_loo)) %>% pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_qc_cov_loo)
    tibw_grp_intensities = tib_input %>% filter(sample_id %in% sid & !is.na(intensity_qc_basic)) %>% pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_qc_basic)
    # if some sample fell out during qc filter, remove from same array
    sid = intersect(sid, colnames(tibw_grp_intensities))
    # skip if <n samples in group. analogous to the generic CoV plot function, but now we require 4 samples per group (after all, we need 3 samples to compute CoV and below we remove 1)
    if(length(sid) < 4) {
      append_log(sprintf("No data available for CoV leave-one-out computation in sample group '%s', skipping plots", grp), type = "info")
      next
    }

    # from tibble to matrix
    mat_grp_intensities = as_matrix_except_first_column(tibw_grp_intensities)
    # # normalize by mode (normalization should have been done upstream already, but only takes a few seconds per group anyway)
    # mat_grp_intensities = normalize_matrix(mat_grp_intensities, group_by_cols = rep(1L, ncol(mat_grp_intensities)), algorithm = "vwmb")

    dropcols = NULL
    mat_loo_cov = matrix(NA, nrow=nrow(tibw_grp_intensities), ncol=length(sid), dimnames = list(tibw_grp_intensities$peptide_id, sid))
    for(sid_exclude in sid) {
      # analogous to the generic CoV plot function, but now we use all samples in group minus the 'leave-one-out'
      # m = as.matrix(tibw_grp_intensities %>% select(!!(setdiff(sid, sid_exclude))))
      m = mat_grp_intensities[ , setdiff(sid, sid_exclude), drop=F]
      rows_fail = matrixStats::rowSums2(!is.na(m)) < 3
      # less than 50 peptides have a value, not a meaningful set of datapoints for CoV analysis
      if(sum(!rows_fail) < 50) {
        dropcols = c(dropcols, sid_exclude)
        next
      }
      m_subset_norm = m[!rows_fail, ]
      # optionally, normalize by mode after removing sample s. This significantly slows down the analysis! For a dataset of N samples, we'd have to normalize the dataset N times (so 200 sample dataset = ~30 minutes for this step alone on a fast workstation and fast vwmb normalization)
      # m_subset_norm = normalize_matrix(m[!rows_fail, ], group_by_cols = rep(1L, ncol(m)), algorithm = "vwmb")

      mat_loo_cov[!rows_fail, sid_exclude] = coefficient_of_variation_vectorized(log2_to_ln(m_subset_norm))
    }
    # debug/QC: bp=boxplot(mat_loo_cov[,order(apply(mat_loo_cov,2,median,na.rm=T))]*100, outline=F, las=2); boxplot_add_text(bp, cex=.5)

    # deal with samples that failed loo computation (eg; too few datapoints after n-value-per-row filter)
    mat_loo_cov = mat_loo_cov[ , !(colnames(mat_loo_cov) %in% dropcols), drop=F]
    if(ncol(mat_loo_cov) == 0) {
      append_log(sprintf("No data available for CoV leave-one-out computation in sample group '%s', skipping plots", grp), type = "info")
      next # skips the bind_rows() below, so no data added to overall data tibble `tib_loo_cov`
    }

    # cov as percentages in all downstream analyses
    mat_loo_cov = mat_loo_cov * 100

    # store median value per sample for downstream plots
    grp_tib_loo = tibble(sample_id = colnames(mat_loo_cov), `median CoV` = matrixStats::colMedians(mat_loo_cov, na.rm = T)) %>%
      arrange(`median CoV`) %>%
      mutate(order_score = dplyr::row_number())

    tib_loo_cov = bind_rows(tib_loo_cov, grp_tib_loo)


    tib_plot = matrix_to_long(mat_loo_cov, value_name = "cov", column_name = "sample_id", row_name = "peptide_id") %>%
      left_join(samples %>% select(sample_id, shortname, group, exclude), by="sample_id")

    tib_plot$exclude = factor(tib_plot$exclude, levels = c(FALSE, TRUE))

    # plot
    p = ggplot(tib_plot, aes(x=cov, colour=shortname, linetype=exclude)) +
      stat_density(geom="line", bw = param_density_bandwidth, adjust = param_density_adjust, position="identity", na.rm = T) +
      scale_linetype_manual(values = c("FALSE"="solid", "TRUE"="dashed")) +
      facet_grid(~group) +
      guides(linetype = "none", #guide_legend(direction = "horizontal", title.position = "top"),
             colour = guide_legend(title = NULL)) +
      # xlim(quantile(tib_plot$cov_log10, probs = c(.001,.999), na.rm=T)) +
      xlim(c(0, quantile(tib_plot$cov, probs = 0.95, na.rm=T))) +
      labs(x="CoV in %, after removing a sample") +
      # labs(title="Effect of removing a sample prior to CoV computation on within-group CoV\nlower value = better CoV after removing sample s", x="log10 Coefficient of Variation") +
      theme_bw() +
      theme(plot.title = element_text(size=10),
            legend.text = element_text(size = ifelse(length(sid) < 4, 10, ifelse(length(sid) < 6, 8, 6)) ),
            legend.position = "bottom", legend.title = element_blank())


    ## split legend into separate plot if more than N entries
    if(dplyr::n_distinct(tib_plot$shortname) > 12) {
      p_split = ggplot_split_legend(p)
      plotlist[[grp]] = p_split$plot
      plotlist[[paste0(grp, "_legend")]] = p_split$legend # + guides(colour = guide_legend(ncol = 2))
    } else {
      plotlist[[grp]] = p
    }


    # debug/QC: sort by CoV, then fix sample_id order by conversion to factor with levels in respective order
    # ggplot(tib_plot %>% arrange(`median CoV`) %>% left_join(samples %>% select(sample_id, shortname, group, exclude), by="sample_id") %>%  mutate(shortname = factor(shortname, levels=shortname)),
    #        aes(x=`median CoV`, y=shortname, color=exclude)) + geom_point() + labs(title=paste("sample group:", grp), subtitle = "Effect on within-group CoV of removing a sample prior to CoV computation", y="") + theme_bw() + theme(legend.position = "bottom")
  }

  if(length(tib_loo_cov) == 0) {
    return(list())
  }


  ### merge data from all 'by group' analyses
  tib_plot = tib_loo_cov %>%
    left_join(samples %>% select(sample_id, shortname, group, exclude), by="sample_id") %>%
    mutate(group = factor(group, levels = rev(unique(samples$group))))

  tib_plot_colors = array(unique(samples_colors$group), dimnames = list(unique(samples$group)))

  # scatterplot
  p_loo_cov = ggplot(tib_plot, aes(x=`median CoV`, y=group, colour=group, shape=ifelse(exclude, "0", "16"))) +
    geom_jitter(height = 0.1, width = 0, alpha = 1) +
    ggrepel::geom_text_repel(data = tib_plot %>% filter(order_score<=2), aes(x=`median CoV`, y=group, colour=group, label=shortname), vjust = .9, hjust = .1, size = 2, max.overlaps = Inf, show.legend = FALSE) +
    scale_shape_manual(values = c("0" = 0, "16" = 16)) +
    guides(shape="none", colour = guide_legend(title = NULL, byrow = T)) +
    # explicitly name all color scale properties to enforce the label order
    scale_color_manual(breaks=names(tib_plot_colors), values = tib_plot_colors, labels=names(tib_plot_colors), aesthetics = c("colour", "fill")) +
    labs(title="Effect of removing a sample prior to CoV computation on within-group CoV\nlower value = better CoV after removing sample s", y="", x="median Coefficient of Variation (%)") +
    theme_bw() +
    theme(plot.title = element_text(size=10), legend.position = "bottom", legend.title = element_blank())

  # summary stats we can print as a table, hard to glance from the figures
  tbl_loo_cov = tib_plot %>%
    select(shortname, group, exclude, `median CoV`) %>%
    mutate(`median CoV` = sprintf("%.1f", `median CoV`))

  append_log_timestamp("leave-one-out CoV plot computations", start_time)
  return(list(loo_bygroup=plotlist, loo_combined=p_loo_cov, tbl_loo_cov=tbl_loo_cov))
}



#' placeholder title
#' @param tib_input todo
#' @param samples todo
#' @param samples_colors todo
#'
#' @importFrom ggpubr theme_pubr
#' @importFrom matrixStats rowSums2
ggplot_coefficient_of_variation = function(tib_input, samples, samples_colors) {

  ## for CoV computation, we need natural log while intensities are log2; log_b(x) = log_d(x) / log_d(b)  @  https://www.purplemath.com/modules/logrules5.htm
  # toy example; given y = log2(x=100), we need z = log10(x);
  # x = 100
  # y = log2(x)
  # z = y / log2(10)
  # x;y;z

  # !! here we use the by-group filtering, and normalization, as configured by the user
  # this is important, because minpep, topN and different normalization make it different from leave-one-out CoV (there, we must use 'mode' norm for speed as we have to normalize the dataset as often as there are samples)
  tibw_abundance_naturallog = tib_input %>%
    select(peptide_id, sample_id, intensity_by_group) %>% # technically, we don't need this, but here to explicitly state input data for now. can this comment out
    filter(!is.na(intensity_by_group)) %>% # this drops all peptides not passing the filter in any sample/group, which makes downstream wide format intensity matrix much smaller
    mutate(intensity_by_group = log2_to_ln(intensity_by_group)) %>%
    pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_by_group)


  ## compute stats per group
  groups = unique(samples$group)
  mat_cov = matrix(NA, nrow=nrow(tibw_abundance_naturallog), ncol=length(groups), dimnames = list(tibw_abundance_naturallog$peptide_id, groups))
  dropcols = NULL
  for(grp in unique(samples$group)) {
    # samples for current group
    sid = intersect(samples %>% filter(group == grp) %>% pull(sample_id),
                    colnames(tibw_abundance_naturallog))
    # skip if <n samples in group
    if(length(sid) < 3) {
      dropcols = c(dropcols, grp)
      append_log(sprintf("no CoV computation for sample group '%s', require at least 3 replicates", grp), type = "info")
      next
    }

    # fast CoV computation
    m = as.matrix(tibw_abundance_naturallog %>% select(!!sid))
    rows_fail = matrixStats::rowSums2(!is.na(m)) < 3
    # less than 50 peptides have a value, not a meaningful set of datapoints for CoV analysis
    if(sum(!rows_fail) < 50) {
      dropcols = c(dropcols, sid_exclude)
      next
    }
    # m[rows_fail, ] = NA # remove rows with less than 3 values (can technically calculate on 2 values, but we chose to require at least 3)
    mat_cov[!rows_fail,grp] = coefficient_of_variation_vectorized(m[!rows_fail,])
  }

  # remove sample groups that have too few replicates
  mat_cov = mat_cov[ , !(colnames(mat_cov) %in% dropcols), drop=F]
  if(ncol(mat_cov) == 0) {
    append_log("No data available for CoV computation, skipping plots", type = "info")
    return(list())
  }

  # cov as percentages in all downstream analyses
  mat_cov = mat_cov * 100
  # debug/QC: bp=boxplot(mat_cov, outline=F, las=2); boxplot_add_text(bp, cex=.5)


  ## summary stats: boxplot. we plot these as text labels onto the ggplot downstream
  bp = boxplot(mat_cov, plot = F)
  cov_data_summ = data.frame(bp$stats)
  colnames(cov_data_summ) = bp$names
  # We split the text labels for the CoV boxplot (all minus last  vs  last) so we can apply different vertical justification
  cov_data_summ_long_14 = cov_data_summ[c(1,4), , drop = F] %>% gather(group, summ) %>% rename(cov = summ)
  cov_data_summ_long_23 = cov_data_summ[c(2,3), , drop = F] %>% gather(group, summ) %>% rename(cov = summ)
  cov_data_summ_long_5 = cov_data_summ[5,,drop = F] %>% gather(group, summ) %>% rename(cov = summ)

  ## summary stats: number of peptides used in CoV computation per group
  tib_cov_group_n = tibble::enframe(colSums(!is.na(mat_cov)), name="group", value = "n")
  # tib_cov_group_n$cov_median = as.numeric(cov_data_summ[3, match(colnames(cov_data_summ), tib_cov_group_n$group)]) # coordinates for plotting at boxplot median line


  # adjust text size to the number of groups
  n_groups = ncol(mat_cov)
  txt_size = 4 - 2*min(10,n_groups) / 10

  # prepare plot data in long format
  tib_plot = matrix_to_long(mat_cov, value_name = "cov", column_name = "group", row_name = "peptide_id") %>%
    mutate(group = factor(group, levels = colnames(mat_cov)))

  ## boxplot
  p_boxplot = ggplot(tib_plot, aes(group, cov, fill = group)) +
    geom_boxplot(outlier.shape = NA, na.rm = T) +
    geom_text(data = cov_data_summ_long_14, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = -.5, nudge_x = .05, size = txt_size) +
    geom_text(data = cov_data_summ_long_23, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0.5, vjust = -.5, size = txt_size) +
    geom_text(data = cov_data_summ_long_5, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = 1, nudge_x = .05, size = txt_size) +
    geom_text(data = tib_cov_group_n, aes(x = group, y = -4, label = n), hjust = 0.5, vjust = 0.5, size = txt_size * 0.8) +
    # geom_text(data = cov_data_summ_long_2, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = .5, nudge_x = .05, size = txt_size) +
    # geom_text(data = tib_cov_group_n, aes(x = group, y = cov_median, label = paste0("n=",n)), hjust = 1, vjust = -.5, nudge_x = -.05, size = txt_size) +
    scale_color_manual(values = array(unique(samples_colors$group), dimnames = list(unique(samples$group))), aesthetics = c("fill")) +
    coord_cartesian(ylim = c(-5, max(cov_data_summ)), clip = 'off') +
    scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(cov_data_summ, na.rm=T) / 10) * 10, by = 10)) +
    labs(x = "", y = "Coefficient of Variation (%)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=ifelse(n_groups<8, 11, 8)))


  ## violin plot
  p_violin = ggplot(tib_plot, aes(group, cov, fill = group)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = .7, na.rm = T) +
    scale_color_manual(values = array(unique(samples_colors$group), dimnames = list(unique(samples$group))), aesthetics = c("fill")) +
    ylim(c(0, max(cov_data_summ))) +
    labs(x = "", y = "Coefficient of Variation (%)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=ifelse(n_groups<8, 11, 8)))


  return(list(violin = p_violin, boxplot = p_boxplot))
}
