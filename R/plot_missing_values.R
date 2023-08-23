
#' placeholder title
#' @param peptides todo
#' @param samples todo
#'
#' @importFrom viridis scale_fill_viridis
#' @importFrom ggpubr theme_pubr
ggplot_peptide_detect_frequency = function(peptides, samples) {
  ## peptide detect counts, mapped to samples
  tib = peptides %>% filter(detect) %>% select(peptide_id, sample_id) %>% add_count(peptide_id, name = "z")

  tib_plot = tib %>%
    count(sample_id, z) %>%
    arrange(desc(z))

  ## sample sorting
  tib_samples = tib_plot %>%
    count(sample_id, wt=n) %>%
    arrange(n) %>%
    select(-n) %>%
    left_join(samples %>% select(sample_id, shortname), by="sample_id")

  tib_plot = tib_plot %>% left_join(tib_samples, by="sample_id")
  tib_plot$sample_id = factor(tib_plot$sample_id, levels = tib_samples$sample_id)
  tib_plot$shortname = factor(tib_plot$shortname, levels = tib_samples$shortname)

  # enforce integer breaks on color scale; https://stackoverflow.com/a/44886993
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]

  p = ggplot(tib_plot, aes(x=shortname, y=n, fill=z)) +
    geom_bar(position="stack", stat="identity", ) +
    coord_flip() +
    scale_fill_distiller(palette = "Spectral", direction=-1, breaks = int_breaks) +
    # viridis::scale_fill_viridis(option = "B") +
    labs(title = "Number of samples in which a peptide is identified vs presence in individual sample", y = "identified peptides", x = "", fill = "#samples") +
    ggpubr::theme_pubr() +
    theme(plot.title = element_text(size=10), legend.position = "right")

  if(nrow(tib_plot) > 50)
    p = p + theme(axis.text.y.left = element_text(size=6))

  return(p)
}



#' placeholder title
#' @param peptides todo
#' @param samples todo
#' @param include_quant todo
#' @param remove_exclude_samples todo
#'
#' @importFrom scales percent_format
ggplot_peptide_detect_frequency_distribution = function(peptides, samples, include_quant = TRUE, remove_exclude_samples = TRUE) {
  # count in how many samples each peptide is identified/quantified, excluding samples flagged by user
  tib = peptides %>% filter(is.finite(intensity))
  if(remove_exclude_samples) {
    tib = tib %>% filter(sample_id %in% (samples %>% filter(!exclude) %>% pull(sample_id)))
  }

  n_samples = n_distinct(tib$sample_id)

  tib = tib %>%
    group_by(peptide_id) %>%
    summarise(identified=sum(detect), quantified=n())

  # summarize cumulative amounts, eg; how many peptides identified in at least N samples ?
  tib_plot = NULL
  for(n in 1:n_samples) {
    tib_plot = bind_rows(tib_plot, tibble(n=n, count=sum(tib$identified >= n), type="identified"))
    tib_plot = bind_rows(tib_plot, tibble(n=n, count=sum(tib$quantified >= n), type="quantified"))
  }

  # besides absolute counts, add the fraction of counts
  tib_plot = tib_plot %>% mutate(n_frac = n/max(n))

  if(!include_quant) {
    tib_plot = tib_plot %>% filter(type == "identified")
  }

  # percentage axis by default has major/minor breaks 0,25,50,75,100%
  n_breaks_perc = .25 * 0:4 # alternatively, additional breaks: .125 * 0:8
  sec_axis_breaks = ceiling(n_breaks_perc * max(tib_plot$n))


  # base ggplot
  p = ggplot(tib_plot, aes(count, n_frac, colour=type))
  if(nrow(tib_plot) > 1) {
    p = p + geom_line()
  }
  p = p +
    geom_point() +
    scale_y_continuous(labels = scales::percent_format(), limits=c(0,1),
                       # sec.axis = dup_axis(name = "Number of samples", labels = function(x) {x * max(tib_plot$n)} )) +
                       sec.axis = sec_axis(~. , breaks = sec_axis_breaks/max(sec_axis_breaks), labels = sec_axis_breaks, name = "Number of samples") ) +
    coord_cartesian(xlim = c(0, max(tib_plot$count))) +
    scale_colour_discrete(labels=c(identified="identified & quantified", quantified="quantified")) +
    labs(title = "data completeness in entire dataset", y = "Fraction of samples", x = "Cumulative amount of peptides") +
    theme_bw() +
    theme(plot.title = element_text(size=10), legend.position = "bottom", legend.title = element_blank())

  # points of interest; at least 50% or 90% of samples
  i5 = head(which(tib_plot$type=="identified" & tib_plot$n_frac>= 0.5), 1)
  i9 = head(which(tib_plot$type=="identified" & tib_plot$n_frac>= 0.9), 1)
  if(length(i5) == 1) {
    p = p + ggrepel::geom_text_repel(aes(label=count), data=tib_plot[i5,], size = 3.5, color="black", direction = "x", segment.alpha = .3, min.segment.length = unit(0.2, 'lines'), max.iter = 10000, seed = 123, box.padding = 0.2) # optionally, only nudge labels in horizontal direction
    # p = p + annotate("text", size = 3.5, label = tib_plot$count[i5], x=tib_plot$count[i5], y = tib_plot$n_frac[i5], hjust=ifelse(tib_plot$count[i5]/max(tib_plot$count) > 0.15, 1.2, 0.5))
    # p = p + annotate("text", label = paste0(tib_plot$count[i5], ", ", tib_plot$n[i5]), x=tib_plot$count[i5], y = tib_plot$n_frac[i5], hjust=1.1)
  }
  if(length(i9) == 1) {
    p = p + ggrepel::geom_text_repel(aes(label=count), data=tib_plot[i9,], size = 3.5, color="black", direction = "x", segment.alpha = .3, min.segment.length = unit(0.2, 'lines'), max.iter = 10000, seed = 123, box.padding = 0.2)
    # p = p + annotate("text", size = 3.5, label = tib_plot$count[i9], x=tib_plot$count[i9], y = tib_plot$n_frac[i9], hjust=ifelse(tib_plot$count[i9]/max(tib_plot$count) > 0.15, 1.2, 0.5))
    # p = p + annotate("text", label = paste0(tib_plot$count[i9], ", ", tib_plot$n[i9]), x=tib_plot$count[i9], y = tib_plot$n_frac[i9], hjust=1.1)
  }

  # analogous for quantified peptide counts
  if(include_quant) {
    q5 = head(which(tib_plot$type=="quantified" & tib_plot$n_frac>= 0.5), 1)
    q9 = head(which(tib_plot$type=="quantified" & tib_plot$n_frac>= 0.9), 1)
    if(length(q5) == 1 && tib_plot$count[i5] != tib_plot$count[q5]) {
      p = p + ggrepel::geom_text_repel(aes(label=count), data=tib_plot[q5,], size = 3.5, color="black", direction = "x", segment.alpha = .3, min.segment.length = unit(0.2, 'lines'), max.iter = 10000, seed = 123, box.padding = 0.2)
      # p = p + annotate("text", size = 3.5, label = tib_plot$count[q5], x=tib_plot$count[q5], y = tib_plot$n_frac[q5], hjust=ifelse(tib_plot$count[q5]/max(tib_plot$count) < 0.85, -0.2, 0.5))
      # p = p + annotate("text", label = paste0(tib_plot$count[q5], ", ", tib_plot$n[q5]), x=tib_plot$count[q5], y = tib_plot$n_frac[q5], hjust=-0.1)
    }
    if(length(q9) == 1 && tib_plot$count[i9] != tib_plot$count[q9]) {
      p = p + ggrepel::geom_text_repel(aes(label=count), data=tib_plot[q9,], size = 3.5, color="black", direction = "x", segment.alpha = .3, min.segment.length = unit(0.2, 'lines'), max.iter = 10000, seed = 123, box.padding = 0.2)
      # p = p + annotate("text", size = 3.5, label = tib_plot$count[q9], x=tib_plot$count[q9], y = tib_plot$n_frac[q9], hjust=ifelse(tib_plot$count[q9]/max(tib_plot$count) < 0.85, -0.2, 0.5))
      # p = p + annotate("text", label = paste0(tib_plot$count[q9], ", ", tib_plot$n[q9]), x=tib_plot$count[q9], y = tib_plot$n_frac[q9], hjust=-0.1)
    }
  }

  return(p)

  ## reference plot code
  # plot(tib_plot$cumsum, tib_plot$n, type="l")
  # j = match(ceiling(0.9 * nrow(samples)), tib_plot$n)
  # text(tib_plot$cumsum[j], tib_plot$n[j], labels = tib_plot$cumsum[j], adj = c(0,0), col=2)
  # j = match(ceiling(0.5 * nrow(samples)), tib_plot$n)
  # text(tib_plot$cumsum[j], tib_plot$n[j], labels = tib_plot$cumsum[j], adj = c(1,1), col=4)
  ## analogous plot in MSnbase package
  # x = msqrobsum::peptide_intensities
  # plotNA(x, pNA = 1/2)
}

# ggplot_peptide_detect_frequency_distribution_v1 = function(peptides, require_detect = TRUE) {
#   tib = peptides %>% filter(!require_detect | detect) %>% count(peptide_id)
#
#   tib_plot = tibble(n = max(tib$n):1, cumsum = 0)
#   for(i in 1:nrow(tib_plot)) {
#     tib_plot$cumsum[i] = sum(tib$n >= tib_plot$n[i])
#   }
#   tib_plot = bind_rows(tibble(n=max(tib_plot$n), cumsum=0), tib_plot) %>% mutate(n_frac = n/max(n))
#
#
#
#   j5 = tail(which(tib_plot$n_frac>= 0.5), 1)
#   j9 = tail(which(tib_plot$n_frac>= 0.9), 1)
#
#   p = ggplot(tib_plot, aes(cumsum, n_frac)) +
#     geom_line() +
#     geom_point() +
#     geom_point(aes(x=tib_plot$cumsum[j5], y = tib_plot$n_frac[j5], color = "red")) +
#     geom_point(aes(x=tib_plot$cumsum[j9], y = tib_plot$n_frac[j9], color = "blue")) +
#     geom_text(aes(x=tib_plot$cumsum[j5], y = tib_plot$n_frac[j5], label = paste0(tib_plot$cumsum[j5], ", ", tib_plot$n[j5]), color = "red"), hjust=ifelse(tib_plot$cumsum[j5] > max(tib_plot$cumsum)/2, 1.1, -0.1)) +
#     geom_text(aes(x=tib_plot$cumsum[j9], y = tib_plot$n_frac[j9], label = paste0(tib_plot$cumsum[j9], ", ", tib_plot$n[j9]), color = "blue"), hjust=ifelse(tib_plot$cumsum[j9] > max(tib_plot$cumsum)/2, 1.1, -0.1)) +
#     scale_y_continuous(labels = scales::percent_format(), limits=c(0,1)) +
#     labs(title = "", y = "Fraction of samples", x = "Cumulative amount of peptides") +
#     theme_bw() +
#     theme(title = element_text(size=8), legend.position = "none")
#   return(p)
# }
