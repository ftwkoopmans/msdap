#' placeholder title
#' @param peptides todo
#' @param samples todo
#' @param isdia todo
#' @importFrom data.table setorder
#' @importFrom patchwork wrap_plots
plot_retention_time_v2 = function(peptides, samples, isdia) {
  param_density_bandwidth = "sj"
  param_density_adjust = 1.0 # alternatively, 0.9

  # TODO: input validation; are all expected columns present in peptides and samples ?
  if(!"intensity_qc_basic" %in% names(peptides)) {
    append_log("plot_retention_time: `intensity_qc_basic` column is missing from peptide tibble", type = "warning")
    return(list())
  }
  start_time = Sys.time()

  # if not DIA, take all peptides. for DIA, only detected (eg; Spectronaut qvalue < x)
  peptides = peptides %>% filter(is.finite(rt) & !is.na(intensity_qc_basic) & (!isdia | detect))

  # samples without replicates cannot be used for this analysis (eg; what is their within-group abundance variation ?)
  key_samples_not_plotted = setdiff(samples$key_sample, peptides$key_sample)
  if(length(key_samples_not_plotted) > 0) {
    append_log(paste("cannot create retention time QC plot for these samples, no peptides with RT and intensity data available (at intensity_qc_basic filter criteria); ", paste(samples$shortname[match(key_samples_not_plotted, samples$key_sample)], collapse = ", ")), type = "warning")
  }

  # never use non-detects to get the reference RT. For DIA, don't use non-detects in these figures
  peptides$temp_rt_nodetect = peptides$rt
  peptides$temp_rt_nodetect[!peptides$detect] = NA
  peptides$temp_int = peptides$intensity_qc_basic

  # use data.tables to summarize the data (reasonably fast, don't need to convert to wide)
  DT = data.table::setDT(peptides)[ , `:=` (temp_rt_diff = rt - stats::median(temp_rt_nodetect, na.rm=T),
                                            temp_int_diff_overall = temp_int - stats::median(temp_int, na.rm=T)), by=key_peptide][ , temp_int_diff := (temp_int - mean(temp_int, na.rm=T)), by=key_peptide_group]

  ## DEBUG: remove data in a few RT bins to simulate spray issues  -->>  how does the plot hold up?
  # DT$rt[DT$key_sample == 2 & DT$rt > 50 & DT$rt <60 ] = NA; print("********* manually discarding data for code debugging")

  ### for testing/debugging; validate raw data
  # hist(peptides %>% filter(key_sample==1 & detect) %>% pull(rt), breaks=80)
  # ggplot(peptides %>% filter(key_sample==1), aes(x = rt, y = temp_rt_diff)) +
  #   geom_point(alpha=0.3, na.rm = T) +
  #   geom_smooth( level=.99, span = 0.01, se=T, na.rm = T) +
  #   # geom_quantile(method = "rqss", lambda = 1, na.rm=T, quantiles=c(.05,.5,.95)) + # slow and non-smoothed
  #   # stat_quantile(formula=y~x, quantiles=c(0.1, 0.5, 0.9), na.rm = T, colour = "red") +  # linear, we don't want this
  #   xlim(overall_rt_xlim) +
  #   coord_cartesian(ylim = c(-2,2))
  # ggplot(tib_sample, aes(x = rt, y = temp_int_diff)) +
  #   geom_point(alpha=0.3, na.rm = T) +
  #   geom_smooth( level=.99, span = 0.01, se=T, na.rm = T) +
  #   xlim(overall_xlim) +
  #   coord_cartesian(ylim = c(-2,2)) # when using ylim() instead of coord_cartesian(), geom_ribbon removes segments beyond the plot limit (which we don't want)



  # bin retention times
  overall_rt_xlim = quantile(DT$rt, probs = c(0.005, 0.995), na.rm = T)
  rtbins = make_bins(DT$rt, from=floor(overall_rt_xlim[1]), to=ceiling(overall_rt_xlim[2]), length = 100)
  rtbins_bin_ids = attr(rtbins, "bin_ids")
  # rtbins_bin_names = attr(rtbins, "bin_names")
  rtbins_bin_means = attr(rtbins, "bin_means")
  DT$temp_rt_bin = as.numeric(rtbins)
  DT$temp_rt_bin_mean = as.numeric(names(rtbins))
  rm(rtbins)

  # count the bin sizes and their respective overall means
  DT_binned = DT[ ,  .(rt_mean = temp_rt_bin_mean[1],
                       size = .N,
                       rt_quantup = quantile(temp_rt_diff, .95, na.rm=T),
                       rt_quantlow = quantile(temp_rt_diff, .05, na.rm=T),
                       rt_median = stats::median(temp_rt_diff, na.rm=T),
                       int_quantup = quantile(temp_int_diff, .95, na.rm=T),
                       int_quantlow = quantile(temp_int_diff, .05, na.rm=T),
                       int_median = stats::median(temp_int_diff, na.rm=T),
                       int_quantup_overall = quantile(temp_int_diff_overall, .95, na.rm=T),
                       int_quantlow_overall = quantile(temp_int_diff_overall, .05, na.rm=T),
                       int_median_overall = stats::median(temp_int_diff_overall, na.rm=T) ),
                  by = .(rt_bin=temp_rt_bin, key_sample)] # note the rename while grouping


  # we created some temp columns by reference, no longer needed as we are working with binned data from here on
  DT[ ,`:=`(temp_rt_nodetect = NULL, temp_int = NULL, temp_rt_diff = NULL, temp_int_diff_overall = NULL, temp_int_diff = NULL, temp_rt_bin = NULL, temp_rt_bin_mean = NULL)]


  # force the creation of empty bins for 'missing bins'. eg; if in some sample at some RT interval no peptides were observed (that pass qc peptide filter), these bins are empty
  list_newdata = list()
  for(key_s in samples$key_sample) {
    DT_binned[key_sample==key_s, ]
    # find empty bins
    empty_bins = setdiff(rtbins_bin_ids, DT_binned[key_sample==key_s, ]$rt_bin)
    # if any, create mock data
    if(length(empty_bins) > 0) {
      list_newdata[[length(list_newdata)+1]] = data.table::data.table(rt_bin = empty_bins, rt_mean = rtbins_bin_means[empty_bins], key_sample=key_s, size=0, rt_quantup=NA, rt_quantlow=NA, rt_median=NA, int_quantup=NA, int_quantlow=NA, int_median=NA, int_quantup_overall=NA, int_quantlow_overall=NA, int_median_overall=NA)
    }
  }


  if(length(list_newdata) > 0) {
    list_newdata[[length(list_newdata)+1]] = DT_binned
    DT_binned = data.table::rbindlist(list_newdata, use.names = T, fill = T)
    rm(list_newdata)
  }




  # reference value for bin sizes
  DT_binned[ , size_overall_median := stats::median(size, na.rm=T), by=rt_bin]

  # smoothed data
  DT_binned[ ,`:=`(rt_quantup_smooth = smooth_loess_custom(rt_bin, rt_quantup, span=.1),
                   rt_quantlow_smooth = smooth_loess_custom(rt_bin, rt_quantlow, span=.1),
                   rt_median_smooth = smooth_loess_custom(rt_bin, rt_median, span=.1),

                   int_quantup_smooth = smooth_loess_custom(rt_bin, int_quantup, span=.1),
                   int_quantlow_smooth = smooth_loess_custom(rt_bin, int_quantlow, span=.1),
                   int_median_smooth = smooth_loess_custom(rt_bin, int_median, span=.1),

                   int_quantup_smooth_overall = smooth_loess_custom(rt_bin, int_quantup_overall, span=.1),
                   int_quantlow_smooth_overall = smooth_loess_custom(rt_bin, int_quantlow_overall, span=.1),
                   int_median_smooth_overall = smooth_loess_custom(rt_bin, int_median_overall, span=.1),

                   size_smooth = smooth_loess_custom(rt_bin, size, span=.1),
                   size_overall_median_smooth = smooth_loess_custom(rt_bin, size_overall_median, span=.1)
                   ),
             by = key_sample]
  # catch smoothing problems for bin size (eg; empty bins)
  DT_binned$size_smooth[!is.finite(DT_binned$size) | DT_binned$size <= 0] = 0

  # sort
  data.table::setorder(DT_binned, rt_mean, na.last = T)

  ### for testing/debugging; plot binned data as-is
  # ggplot(as_tibble(DT_binned[key_sample==1, ]), aes(x = rt_mean, y = size)) +
  #   geom_line(na.rm = T, colour = "black") +
  #   geom_line(aes(y=size_overall_median), na.rm = T, colour = "grey")
  # ggplot(as_tibble(DT_binned[key_sample==1, ]), aes(x = rt_mean, y = rt_median)) +
  #   geom_line(na.rm = T, colour = "black") +
  #   geom_line(aes(y=rt_quantup), na.rm = T, colour = "blue") +
  #   geom_line(aes(y=rt_quantlow), na.rm = T, colour = "blue") +
  #   coord_cartesian(ylim = c(-1,1))
  # ggplot(as_tibble(DT_binned[key_sample==1, ]), aes(x = rt_mean, y = int_median)) +
  #   geom_line(na.rm = T, colour = "black") +
  #   geom_line(aes(y=int_quantup), na.rm = T, colour = "blue") +
  #   geom_line(aes(y=int_quantlow), na.rm = T, colour = "blue") +
  #   coord_cartesian(ylim = c(-2,2)) # when using ylim() instead of coord_cartesian(), geom_ribbon removes segments beyond the plot limit (which we don't want)



  append_log_timestamp("RT plots: preparing data", start_time)
  start_time = Sys.time()

  ## all samples in a single plot
  # note; using match for speed
  i = match(peptides$key_sample, samples$key_sample)
  peptides$exclude = samples$exclude[i]
  # update; we now store group as a factor with levels in the same order as input sample table, so plot legends are sorted in 'standard order' as expected from sample metadata table
  peptides$group = droplevels(factor(samples$group[i], levels = unique(samples$group))) #previously; peptides$group = samples$group[i]
  ngroup = nlevels(samples$group) #previously; ngroup = dplyr::n_distinct(samples$group)


  p_all_rt_distributions = ggplot(peptides, aes(x=temp_rt_nodetect, group=sample_id, linetype = exclude)) +
    geom_line(stat="density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm=T, alpha=0.3) +
    xlim(overall_rt_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
    labs(x="Retention time (min)", y="(detected) peptide density") +
    facet_grid(~"all samples, no color-coding") +
    theme_bw() +
    theme(legend.position = "none", legend.title = element_blank())

  p_all_rt_distributions_colour_groups = ggplot(peptides, aes(x=temp_rt_nodetect, group=sample_id, colour=group, linetype = exclude)) +
    geom_line(stat="density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm=T) +
    guides(linetype = "none") +
    xlim(overall_rt_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
    labs(x="Retention time (min)", y="(detected) peptide density") +
    facet_grid(~"all samples, color-coded by group") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = ifelse(ngroup < 4, 10, ifelse(ngroup < 6, 8, 6)) ),
          legend.title = element_text(size=10) )

  p_all_rt_distributions_collapse_groups = ggplot(peptides, aes(x=temp_rt_nodetect, colour=group)) +
    geom_line(stat="density", bw = param_density_bandwidth, adjust = param_density_adjust, na.rm=T) +
    xlim(overall_rt_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
    labs(x="Retention time (min)", y="(detected) peptide density") +
    facet_grid(~"overall density by sample group") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = ifelse(ngroup < 4, 10, ifelse(ngroup < 6, 8, 6)) ),
          legend.title = element_text(size=10) )


  ## color-code by group, each line/sample in different shades of the same color
  # clr_sample_by_group = samples %>% select(sample_id, group) %>% add_column(clr = NA)
  # grps = table(clr_sample_by_group$group)
  # basecolors = colorRampPalette(RColorBrewer::brewer.pal(min(9, length(ugrp)), "Set1"), space = "Lab")(length(ugrp))
  # # clr_ref = samples_colors_long %>% filter(prop=="group") %>% distinct(val, .keep_all = T)
  # # basecolors = clr_ref$clr[match(names(grps), clr_ref$val)]
  # for(i in seq_along(grps)) {
  #   g = names(grps[i])
  #   n = as.numeric(grps[i])
  #   clr = basecolors[i]
  #   if(n > 1) {
  #     clr = colorspace::darken(basecolors[i], amount = seq(from=-0.3, to=.3, length.out=n))
  #   }
  #   clr_sample_by_group$clr[clr_sample_by_group$group==g] = clr
  # }
  # m = nrow(clr_sample_by_group)
  # clr_alpha = ifelse(m<10, "DD", ifelse(m<25, "AA", "99"))
  # clr_sample_by_group$clr = paste0(clr_sample_by_group$clr, clr_alpha)
  #
  # p_all_rt_distributions_colour_groups = p_all_rt_distributions +
  #   scale_colour_manual(values=array(clr_sample_by_group$clr, dimnames = list(clr_sample_by_group$sample_id)))
  # p_all_rt_distributions_colour_groups

  # consistent plot limits
  plotlim_rt_diff = min(5, max(abs(c(quantile(DT_binned$rt_quantlow_smooth, probs = .005, na.rm = T), quantile(DT_binned$rt_quantup_smooth, probs = .995, na.rm = T)))))
  plotlim_rt_diff = c(-1,1) * max(0.1, plotlim_rt_diff) # guard against zero diff
  plotlim_int_diff = c(-1,1) * (min(5, max(abs(c(quantile(DT_binned$int_quantlow_smooth_overall, probs = .005, na.rm = T), quantile(DT_binned$int_quantup_smooth_overall, probs = .995, na.rm = T))))))


  plotlist = list()
  for(key_s in samples$key_sample) { # key_s=2
    s = samples %>% filter(key_sample == key_s)
    DT_binned_sample_s = DT_binned[key_sample==key_s, ]
    DT_binned_sample_s$label = sprintf("%s @ %s%s", s$shortname, s$group, ifelse(s$exclude,"  >>  EXCLUDED  <<",""))

    # p = ggplot(size=1) +
    #   facet_grid(plotlevel~label, scales = "free_y", labeller = labeller(plotlevel=function(x) { c("peptide\ncount", "retention time -\noverall median","log2 intensity -\noverall median","log2 intensity -\ngroup mean")[as.integer(x)] })) +
    #   # geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=1), aes(x=rt_mean, y=size), size=1, na.rm = T, colour="cyan") + # for testing, plot without smoothing
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=1), na.rm=T, aes(x=rt_mean, y=size_overall_median_smooth), size=1, colour="grey", alpha = 0.7) +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=1), na.rm=T, aes(x=rt_mean, y=size_smooth), size=1, colour="#2c7fb8") +
    #   # geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=1), aes(x=rt_mean, y=size_smooth), size=1, na.rm = T, colour="#2c7fb8") +
    #   # geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=1), aes(x=rt_mean, y=size_overall_median_smooth), size=1, na.rm = T, colour="#bbbbbb") +
    #
    #   geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean), y=0, colour="#555555", linetype = "dashed") +
    #   geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean), y=0, colour="#555555", linetype = "dashed") +
    #   geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(x = rt_mean), y=0, colour="#555555", linetype = "dashed") +
    #
    #   # geom_ribbon(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean, ymin = rt_quantlow_smooth, ymax = rt_quantup_smooth), alpha = 0.6, fill = "grey") +
    #   # geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean, y=rt_median_smooth), size=1, colour = "#c51b8a") + # 0.5 + 2*(size_smooth/max(size_smooth))
    #   ## 'stepped area fill' implementation is inspired by this example code; https://stackoverflow.com/a/43621741
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean, y=rt_quantlow_smooth), size=0.5, colour = "grey") +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean, y=rt_quantup_smooth), size=0.5, colour = "grey") +
    #   geom_rect(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(xmin = rt_mean, xmax = lead(rt_mean), ymin = rt_quantlow_smooth, ymax = rt_quantup_smooth), fill = "grey", alpha = 0.4) +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=2), na.rm=T, aes(x = rt_mean, y=rt_median_smooth, size=size_smooth), colour = "black") + # 0.5 + 2*(size_smooth/max(size_smooth))
    #
    #   # geom_ribbon(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean, ymin = int_quantlow_smooth_overall, ymax = int_quantup_smooth_overall), alpha = 0.6, fill = "grey") +
    #   # geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean, y=int_median_smooth_overall), size=1, colour = "#e34a33") +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean, y=int_quantlow_smooth_overall), size=0.5, colour = "grey") +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean, y=int_quantup_smooth_overall), size=0.5, colour = "grey") +
    #   geom_rect(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(xmin = rt_mean, xmax = lead(rt_mean), ymin = int_quantlow_smooth_overall, ymax = int_quantup_smooth_overall), fill = "grey", alpha = 0.4) +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=3), na.rm=T, aes(x = rt_mean, y=int_median_smooth_overall, size=size_smooth), colour = "#e34a33") +
    #
    #   # geom_ribbon(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(x = rt_mean, ymin = int_quantlow_smooth, ymax = int_quantup_smooth), alpha = 0.6, fill = "grey") +
    #   # geom_line(data = DT_binned_sample_s %>% add_column(plotlevel=3), aes(x = rt_mean, y=int_median_smooth), size=1, colour = "#e34a33") +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(x = rt_mean, y=int_quantlow_smooth), size=0.5, colour = "grey") +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(x = rt_mean, y=int_quantup_smooth), size=0.5, colour = "grey") +
    #   geom_rect(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(xmin = rt_mean, xmax = lead(rt_mean), ymin = int_quantlow_smooth, ymax = int_quantup_smooth), fill = "grey", alpha = 0.4) +
    #   geom_step(data = DT_binned_sample_s %>% add_column(plotlevel=4), na.rm=T, aes(x = rt_mean, y=int_median_smooth, size=size_smooth), colour = "#c51b8a") +
    #
    #   scale_size(range = c(0.25, 1.5)) +
    #   guides(size = "none") +
    #   labs(x="Retention time (min)", y="") +
    #   theme_bw()
    # # highlight 'exclude' samples
    # if(s$exclude) {
    #   p = p + theme(strip.background.x = element_rect(fill="#fdbb84"))
    # }


    p1 = ggplot(DT_binned_sample_s) +
      geom_step(na.rm=T, aes(x=rt_mean, y=size_overall_median_smooth), size=1.5, colour="grey", alpha = 0.7) +
      geom_step(na.rm=T, aes(x=rt_mean, y=size_smooth), size=1, colour="black") +
      facet_grid(~ label) +
      scale_size(range = c(0.25, 1.5)) +
      labs(x="", y="peptide count") +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

    p2 = ggplot(DT_binned_sample_s) +
      geom_line(na.rm=T, aes(x = rt_mean), y=0, colour="#555555", linetype = "dashed") +
      geom_step(na.rm=T, aes(x = rt_mean, y=rt_quantlow_smooth), size=0.5, colour = "grey") +
      geom_step(na.rm=T, aes(x = rt_mean, y=rt_quantup_smooth), size=0.5, colour = "grey") +
      geom_rect(na.rm=T, aes(xmin = rt_mean, xmax = lead(rt_mean), ymin = rt_quantlow_smooth, ymax = rt_quantup_smooth), fill = "grey", alpha = 0.4) +
      geom_step(na.rm=T, aes(x = rt_mean, y=rt_median_smooth, size=size_smooth), colour = "#2c7fb8") +
      coord_cartesian(ylim = plotlim_rt_diff) +
      scale_size(range = c(0.25, 1.5)) +
      guides(size = "none") +
      labs(x="", y="retention time -\noverall median") +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

    p3 = ggplot(DT_binned_sample_s) +
      geom_line(na.rm=T, aes(x = rt_mean), y=0, colour="#555555", linetype = "dashed") +
      geom_step(na.rm=T, aes(x = rt_mean, y=int_quantlow_smooth_overall), size=0.5, colour = "grey") +
      geom_step(na.rm=T, aes(x = rt_mean, y=int_quantup_smooth_overall), size=0.5, colour = "grey") +
      geom_rect(na.rm=T, aes(xmin = rt_mean, xmax = lead(rt_mean), ymin = int_quantlow_smooth_overall, ymax = int_quantup_smooth_overall), fill = "grey", alpha = 0.4) +
      geom_step(na.rm=T, aes(x = rt_mean, y=int_median_smooth_overall, size=size_smooth), colour = "#e34a33") +
      coord_cartesian(ylim = plotlim_int_diff) +
      scale_size(range = c(0.25, 1.5)) +
      guides(size = "none") +
      labs(x="Retention time (min)", y="log2 intensity -\noverall median") +
      theme_bw()

    if(s$exclude) {
      p1 = p1 + theme(strip.background.x = element_rect(fill="#fdbb84"))
    }

    plotlist[[length(plotlist) + 1]] = patchwork::wrap_plots(p1, p2, p3, ncol = 1)
    # print(plotlist[[length(plotlist)]])
  }

  append_log_timestamp("RT plots: creating plots", start_time)
  return(list(key_samples_not_plotted = key_samples_not_plotted, rt_distributions_all = p_all_rt_distributions, rt_distributions_colour_groups = p_all_rt_distributions_colour_groups, rt_distributions_collapse_groups = p_all_rt_distributions_collapse_groups, rt_by_sample = plotlist))
}












# ### result 1: plot the binned data as-is
# ggplot(as_tibble(y), aes(x = temp_rt_bin, y=temp_rt_bin_qmed, ymin = temp_rt_bin_qlow, ymax = temp_rt_bin_qup)) +
#   geom_ribbon(alpha = 0.6, fill = "grey") +
#   geom_line(na.rm = T, colour = "black") +
#   theme_light() +
#   coord_cartesian(ylim = c(-2,2)) # when using ylim() instead of coord_cartesian(), geom_ribbon removes segments beyond the plot limit (which we don't want)
#
# ### result 2: manually smooth binned data
# smooth_loess_custom = function(x, y, span=.1) {
#   predict(loess(y ~ x, span = span), x)
# }
#
# p = ggplot(as_tibble(y) %>% mutate(y=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qmed, span=.1),
#                                    ymin=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qlow, span=.1),
#                                    ymax=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qup, span=.1)),
#            aes(x = temp_rt_bin, y=y, ymin = ymin, ymax = ymax)) +
#   geom_ribbon(alpha = 0.6, fill = "grey") +
#   geom_line(na.rm = T, colour = "black") +
#   theme_light() +
#   coord_cartesian(ylim = c(-2,2)) # when using ylim() instead of coord_cartesian(), geom_ribbon removes segments beyond the plot limit (which we don't want)
# p
#
# ### TODO: desireable histogram-line/density to show amount of peptides and grey line for median amount over all samples
# # this should work for samples that have 0 peptides in a bin
# p2 = ggplot(y, aes(x=temp_rt_bin, y=temp_rt_bin_size)) +
#   # geom_step() + # have to smooth due to binning;
#   geom_smooth(span = 0.25, se=T, na.rm = T) +
#   theme_light()
#
# p
# p2
#
#
# ### combine these using facets ?
#
# tib_plot = as_tibble(y) %>% mutate(temp_rt_bin_qmed=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qmed, span=.1),
#                                    temp_rt_bin_qlow=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qlow, span=.1),
#                                    temp_rt_bin_qup=smooth_loess_custom(temp_rt_bin, temp_rt_bin_qup, span=.1),
#
#                                    temp_rt_bin_int_qmed=smooth_loess_custom(temp_rt_bin, temp_rt_bin_int_qmed, span=.1),
#                                    temp_rt_bin_int_qlow=smooth_loess_custom(temp_rt_bin, temp_rt_bin_int_qlow, span=.1),
#                                    temp_rt_bin_int_qup=smooth_loess_custom(temp_rt_bin, temp_rt_bin_int_qup, span=.1),
#                                    temp_rt_bin_refrt = rtbins[temp_rt_bin]
#                                    )
# tib_plot = bind_rows(tib_plot %>% add_column(plotlevel=1),
#                      tib_plot %>% add_column(plotlevel=2),
#                      tib_plot %>% add_column(plotlevel=3))
#
#
# ggplot() +
#   facet_wrap(.~plotlevel, scales = "free_y", dir = "v", labeller = labeller(plotlevel=function(x) { c("number of peptides", "RT - overall median","log2 intensity - group average")[as.integer(x)] })) +
#   geom_smooth(data = subset(tib_plot, plotlevel == 1), aes(x=temp_rt_bin_mean, y=temp_rt_bin_size), size=1, span = 0.2, se=F, na.rm = T, colour="#2ca25f") +
#   geom_line(data = subset(tib_plot, plotlevel != 1), aes(x = temp_rt_bin_mean), y=0, colour="#555555", linetype = "dashed") +
#   geom_ribbon(data = subset(tib_plot, plotlevel == 2), aes(x = temp_rt_bin_mean, ymin = temp_rt_bin_qlow, ymax = temp_rt_bin_qup), na.rm = T, alpha = 0.6, fill = "grey") +
#   geom_line(data = subset(tib_plot, plotlevel == 2), aes(x = temp_rt_bin_mean, y=temp_rt_bin_qmed), size=1, na.rm = T, colour = "#c51b8a") +
#   geom_ribbon(data = subset(tib_plot, plotlevel == 3), aes(x = temp_rt_bin_mean, ymin = temp_rt_bin_int_qlow, ymax = temp_rt_bin_int_qup), na.rm = T, alpha = 0.6, fill = "grey") +
#   geom_line(data = subset(tib_plot, plotlevel == 3), aes(x = temp_rt_bin_mean, y=temp_rt_bin_int_qmed), size=1, na.rm = T, colour = "#e34a33") +
#   labs(x="Retention time", y="") +
#   theme_bw()








############################### plot version one; density scatterplots
#' #' placeholder title
#' #' @param tib_input todo
#' #' @param samples todo
#' #'
#' #' @importFrom ggpointdensity geom_pointdensity
#' #' @importFrom viridis scale_color_viridis
#' plot_retention_time = function(tib_input, samples) {
#'   tib = inner_join(tib_input %>% filter(is.finite(rt) & !is.na(intensity_qc_basic)) %>% select(sample_id, peptide_id, detect, rt, intensity = intensity_qc_basic),
#'                    samples %>% select(sample_id, shortname, group, exclude),
#'                    by = "sample_id")
#'   tib_detect = tib %>% filter(detect)
#'   overall_xlim = quantile(tib$rt, probs = c(0.005, 0.995), na.rm = T)
#'
#'   ## all samples in a single plot
#'   p_all_rt_distributions = ggplot(tib_detect, aes(x=rt, colour=shortname, labels=shortname, linetype = exclude)) +
#'     geom_line(stat="density", na.rm=T) +
#'     xlim(overall_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
#'     labs(x="retention time (min)", y="(detected) peptide density") +
#'     theme_bw() +
#'     theme(legend.position = "none", legend.title = element_blank())
#'
#'   # optionally, scale retention time density by peptide intensity
#'   # tib = tib %>% group_by(sample_id) %>% mutate(intensity_sample_fraction = 2^intensity / sum(2^intensity))
#'
#'   # overall retention-time distributions, per group so plots are readable for large datasets
#'   l_group_rt_distributions = list()
#'   for (g in unique(tib$group)) {
#'     n_samples = sum(samples$group == g)
#'     # weight=intensity_sample_fraction, # optionally, scale retention time density by peptide intensity. https://github.com/tidyverse/ggplot2/issues/2900
#'     p = ggplot(tib_detect %>% filter(group == g), aes(x = rt, colour = shortname, labels = shortname, linetype = exclude)) +
#'       # stat_density(adjust = 1/3, geom="line", position="identity", na.rm=T) + # optionally, adjust the bandwidth/sensitivity
#'       geom_line(stat = "density", na.rm=T) +
#'       xlim(overall_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
#'       labs(x="retention time (min)", y="(detected) peptide density", colour="sample", linetype = "exclude sample") +
#'       guides(linetype = "none", #guide_legend(direction = "horizontal", title.position = "top"),
#'              colour = guide_legend(title = NULL)) +
#'       facet_grid(~group) +
#'       theme_bw() +
#'       theme(legend.position = "bottom",
#'             legend.text = element_text(size = ifelse(n_samples < 4, 10, ifelse(n_samples < 6, 8, 6)) ),
#'             legend.title = element_text(size=10) )
#'
#'     l_group_rt_distributions[[g]] = p
#'   }
#'
#'   # overall reference retention time; only use detected
#'   pep_reference_rt = tib_detect %>% group_by(peptide_id) %>% add_tally() %>% filter(n > 1) %>% summarise(rt_avg = mean(rt, na.rm = T))
#'   # group count -->> require at least 2 samples with quantification -->> store mean intensity per peptide*group
#'   pep_reference_intensity = tib %>% group_by(peptide_id, group) %>% add_tally() %>% filter(n > 1) %>% summarise(int_avg = mean(intensity, na.rm = T))
#'   #
#'   tib_merged = tib %>%
#'     left_join(pep_reference_rt, by = "peptide_id") %>%
#'     left_join(pep_reference_intensity, by = c("peptide_id", "group"))
#'   # scores that we want to show for each peptide
#'   tib_merged = tib_merged %>% filter(is.finite(rt_avg))
#'   tib_merged = tib_merged %>%
#'     add_column(rt_diff = tib_merged$rt - tib_merged$rt_avg) %>%
#'     add_column(int_diff = tib_merged$intensity - tib_merged$int_avg)
#'
#'   plotlist = list()
#'   for (sample_name in samples$shortname) { # sample_name = samples$shortname[1]
#'     tib_sample = tib_merged %>% filter(shortname == sample_name)
#'
#'     # data for current sample:
#'     # RT  *  (RT - `reference RT`)  -->> reference is median in entire dataset
#'     # intensity  *  log2(intensity_norm - `reference intensity_norm`)  -->> reference is within-group mean, only for values with 2+ detect in group
#'     tib_sample_long = bind_rows(
#'       tib_sample %>% select(shortname, group, exclude, rt, detect, value = rt_diff) %>% add_column(key = paste0(sample_name, ": RT - median")),
#'       tib_sample %>% select(shortname, group, exclude, rt, detect, value = int_diff) %>% drop_na() %>% add_column(key = paste0(sample_name, ": log2(int / group mean)"))
#'     )
#'     tib_sample_long$key = factor(tib_sample_long$key, levels = unique(tib_sample_long$key)) # factor with levels in order as we just added them to tib_plot, forcing order for downstream plot
#'
#'     # scatterplot with a facet on the `key`
#'     p = ggplot(tib_sample_long, aes(x = rt, y = value, labels = key)) +
#'       ggpointdensity::geom_pointdensity(adjust = 4, size = .25, alpha = 0.8, na.rm = T) + # aes(shape = detect),
#'       # scale_colour_gradient(low="grey", high = "purple") +
#'       viridis::scale_color_viridis(option = "D", direction = -1) +
#'       xlim(overall_xlim) + # use xlim instead of coord_cartesian to compute the densities only on the subset of data points within this limited RT window, to prevent influence from far outliers
#'       coord_cartesian(ylim = c(-5, 5)) +
#'       scale_shape_manual(values = c(16, 15)) +
#'       scale_y_continuous("", sec.axis = dup_axis(name = ifelse(tib_sample_long$exclude, "EXCLUDED SAMPLE", paste("group:", tib_sample_long$group)))) + # add an axis on the right-hand side with some custom label
#'       labs(x = "retention time (min)") +
#'       guides(colour = "none", shape = guide_legend(title = "detected?")) +
#'       facet_grid(~key) +
#'       # facet_wrap(. ~ key, scales="free_y") +
#'       theme_bw() +
#'       theme(legend.position = "none",
#'             axis.ticks.y.right = element_blank(),
#'             axis.text.y.right = element_blank(),
#'             axis.title.y = element_text(size=8))
#'
#'
#'     ## v1: plot all data points with some alpha. downside; will be very slow for huge datasets with 40k peptides * 2 plots * N samples
#'     # ggplot(tib_sample_long, aes(x = rt, y = value, labels = key)) +
#'     #   geom_point(aes(colour = key, shape = detect), shape=16, alpha = .2) +
#'     #   stat_smooth(method = "gam", formula = y ~ s(x), se = F, colour = "black") +
#'     #   scale_colour_manual(values = c("#1f78b444", "#984ea344")) +
#'     #   coord_cartesian(xlim=overall_xlim, ylim = c(-5, 5)) +
#'     #   scale_shape_manual(values = c(16, 15)) +
#'     #   scale_y_continuous("", sec.axis = sec_axis(~., labels = NULL, breaks = NULL,
#'     #                                              name = ifelse(tib_sample_long$exclude, "EXCLUDED SAMPLE", paste("group:", tib_sample_long$group)))) + # trick to add an axis on the right-hand side with some custom label
#'     #   labs(x = "retention time (min)") +
#'     #   guides(colour = "none", shape = guide_legend(title = "detected?")) +
#'     #   facet_grid(~key) +
#'     #   theme_bw() +
#'     #   theme(legend.position = "none")
#'
#'     plotlist[[sample_name]] = p
#'   }
#'
#'   return(list(rt_distributions_all = p_all_rt_distributions, rt_distributions_bygroup = l_group_rt_distributions, rt_by_sample = plotlist))
#' }
