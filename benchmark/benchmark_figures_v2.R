# plot goal: overall score for normalization
# in each dataset*contrast*algo_de; score per algo_norm = pAUC's % difference from best pAUC (within this ROC)
# x = pAUC performance gap to best-in-dataset
# y = algo_norm
# main figure; all algo_de collapsed per dataset*contrast (use algo_de as shape)
# supplementary; y-axis = dataset*contrast*algo_de
#
#
#
## plot goal: overall score for algo_de
# pAUC to find 'overall best performer' (and conclude that ebayes is slightly behind and msempire and msqrob are similar, with few tossups depending on the dataset)
# 1 normalization algo; vsn&modebetween (or whatever was best overall)
# x = pAUC
# y = dataset*contrast
# color = algo_de
# pch = algo_de
# ggplot(tib_plot, aes(pauc, dataset_contrast, colour=algo_de)) + geom_jitter()
#
#
#
## plot goal: qvalue estimates vary between methods and datasets, as seen in ROC
# number of FPs at some pvalue threshold is quite different (while ROCs are comparable)
# panel A: illustrate how ROC's are all pretty good, but estimates of 'optimal threshold' are not always at qval001
# panel A: zoomed ROC for LFQbench/vsn&mb --> youden index, p-value threshold, number of true/false-positives as text in legenda
#
# panel B: summary statistics when applying this to all datasets
# for each algo_de: FPR at qval001, FPR at youden index, pAUC at youden index
# x = Qvalue threshold at 10% FDR
# y = dataset*contrast
# color = algo_de
# pch = algo_de
## ?? is ROC the same with qvalue as pvalue? otherwise, threshold is the uncorrected pvalue
#
# supplementary: ROC per dataset (2 panels, entire ROC + zoom) at 1 normalization; compare algo_de -->> observe that youden index and qval001 cutoff are sometimes far apart for msqrob/msempire
##
# -->> use youden index in ROC to show eBayes is too conservative and msempire and msqrob, for some datasets, miss the mark in estimating qvalues
# x = algo_de
# y = pvalue at youden index
# data = values from all datasets and contrasts
# color = dataset
# pch = contrast number within dataset (?)


##############################################################################################################################
######################################################## combine data ########################################################
##############################################################################################################################

# rm(list = ls(all.names = TRUE)) # clear everything from memory
# cat("\014") # clear terminal (send the control+L character)
# devtools::load_all() # load our R package
library(tidyverse)
library(pROC)
load("C:/temp/benchmark_datasets.RData")
# load("C:/temp/benchmark_datasets_topN-10.RData")
# load("C:/temp/benchmark_datasets_topN-5.RData")

roc_data = roc_data %>% mutate(algo_lbl = paste0(algo_de, " @ ", algo_norm),
                               dataset_name = ds_name,
                               dataset_contrast_id = paste(ds_name, contrast),
                               dataset_contrast_norm_id = paste(dataset_contrast_id, algo_norm),
                               dataset_contrast_norm_algode_id = paste(dataset_contrast_id, algo_lbl) )


# create pROC object and extract 'youden index' (threshold/TPR/FPR), pAUC at 5% and 10% FDR, FPR/TPR at qvalue<=0.01/0.05
# store pROC object and summary stats in a tibble (long format)
dcna_id = unique(roc_data$dataset_contrast_norm_algode_id)
tib_roc = tibble()
for(id in dcna_id) { #id = grep("lfqbench_.*vsn ebayes", dcna_id, value=T)
  # tib = roc_data %>% filter(is.finite(pvalue) & dataset_contrast_norm_algode_id == id & peptides_used_for_dea >= 10) # TOPN TEST
  tib = roc_data %>% filter(is.finite(pvalue) & dataset_contrast_norm_algode_id == id)

  # if using qvalues and minlog10 values are needed; use the precomputed qvalue_minlog10 as it correctly handles stats methods that yield qvalue=0 (eg; msempire)
  # compute ROC
  roc_obj_pauc = pROC::roc(tib$classification, tib$qvalue, levels=c("background", "foreground"), direction=">", partial.auc=c(1, 0.95), partial.auc.correct=F, partial.auc.focus="specificity")
  pauc = roc_obj_pauc$auc
  roc_obj = pROC::roc(tib$classification, tib$qvalue, levels=c("background", "foreground"), direction=">")
  coord_best = pROC::coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity"), transpose = FALSE, best.method="youden")
  coord_signif = pROC::coords(roc_obj, 0.01, input="threshold", ret=c("threshold", "specificity", "sensitivity", "tp", "fp", "tpr", "fpr"), transpose = FALSE)
  # cannot use to get sensitivity at some specificity/fpr, because it'll yield NA for threshold in some cases (maybe a bug?). Instead, we get closest value from entire list of coordinates
  sc = pROC::coords(roc_obj, input="specificity", ret=c("threshold", "specificity", "sensitivity", "tp", "fp", "tpr", "fpr"), transpose = FALSE)
  coord_fpr5 = sc[tail(which(sc$specificity <= 0.95),1), ]
  coord_fpr10 = sc[tail(which(sc$specificity <= 0.9),1), ]

  names(coord_best) = paste0("best_", names(coord_best))
  names(coord_signif) = paste0("signif_", names(coord_signif))
  names(coord_fpr5) = paste0("fpr5_", names(coord_fpr5))
  names(coord_fpr10) = paste0("fpr10_", names(coord_fpr10))

  tib_roc = bind_rows(tib_roc,
                      as_tibble(c(coord_best, coord_signif, coord_fpr5, coord_fpr10)) %>% add_column(dataset_contrast_norm_algode_id=id, pauc=as.numeric(pauc), proc_pauc = list(roc_obj_pauc), proc = list(roc_obj), .before = 1))
}
# add dataset metadata to ROC tibble
tib_roc = tib_roc %>% inner_join(roc_data %>% select(dataset_contrast_id, algo_lbl, dataset_contrast_norm_algode_id, algo_de, algo_norm, dataset_name, contrast) %>% distinct(dataset_contrast_norm_algode_id, .keep_all=T), by="dataset_contrast_norm_algode_id")




# #### TOPN TEST
# pdf("C:/temp/roc_topn.pdf")
# #### regular plots; compare normalization within same algo_de
# for(id in unique(roc_data$dataset_contrast_id)) {
#   x = tib_roc %>% filter(dataset_contrast_id == id) %>% arrange(desc(pauc))
#   x_clr = rainbow(nrow(x), alpha = 0.8)
#   for(i in 1:nrow(x)) {
#     pROC::plot.roc(x$proc_pauc[[i]], add=(i>1), col=x_clr[i], main=id, lwd=1, cex.main=0.6)
#   }
#   legend("bottomright", legend = sprintf("pAUC:%.2f%% TP:%s FP:%s TPR:%.2g FPR:%.2g %s", x$pauc*100, x$signif_tp, x$signif_fp, x$signif_tpr, x$signif_fpr, x$algo_lbl), col=x_clr, bty='n', cex=.6, lty=1)
#
#   for(alg in unique(roc_data$algo_de)) { #id=roc_data$dataset_contrast_id[1]; alg="ebayes"
#     x = tib_roc %>% filter(dataset_contrast_id == id & algo_de == alg) %>% arrange(desc(pauc))
#     x_clr = rainbow(nrow(x), alpha = 0.8)
#     for(i in 1:nrow(x)) {
#       pROC::plot.roc(x$proc_pauc[[i]], add=(i>1), col=x_clr[i], main=paste(id, alg), lwd=1, cex.main=0.6)
#     }
#     legend("bottomright", legend = sprintf("pAUC:%.2f%% TP:%s FP:%s TPR:%.2g FPR:%.2g %s", x$pauc*100, x$signif_tp, x$signif_fp, x$signif_tpr, x$signif_fpr, x$algo_norm), col=x_clr, bty='n', cex=.6, lty=1)
#   }
# }
#
# dev.off()
#
# # l = x$proc_pauc; names(l) = x$algo_lbl
# # pROC::ggroc(l, alpha=0.8) +
# #   geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed") +
# #   scale_colour_hue() +
# #   labs(x="FPR", y="TPR") +
# #   theme_bw()





##############################################################################################################################
###########################################################  plot  ###########################################################
##############################################################################################################################

pdf("C:/temp/new.pdf")

#### regular plots; foldchange distributions
for(id in unique(roc_data$dataset_contrast_id)) { # id=unique(roc_data$dataset_contrast_id)[3]
  x = roc_data %>% filter(dataset_contrast_id == id & algo_de == "ebayes")
  x$classification = droplevels(factor(x$classification, levels = c("discard", "background", "foreground")))
  x$algo_norm = factor(x$algo_norm, levels = sort(unique(x$algo_norm)))
  x$box = paste(x$classification, x$algo_norm)

  # fix order
  x$box = factor(x$box, levels = as.vector(outer(levels(x$classification), levels(x$algo_norm), paste)) )

  # compute boxplot stats for each, then find the outer limits
  x_limits = x %>% group_by(box) %>% summarise(stats = list(range(boxplot.stats(foldchange.log2)$stats)))
  foldchange_limits = range(unlist(x_limits$stats, recursive = T))

  print(ggplot(x, aes(box, foldchange.log2, fill=classification)) +
          geom_boxplot(outlier.shape = NA) +
          scale_x_discrete(labels=sub(".* ", "", levels(x$box))) +
          labs(y="log2 foldchange between groups in contrast", x="", title = id) +
          coord_flip(ylim = foldchange_limits) +
          theme_bw() + theme(plot.title = element_text(size=8)))
  # print(ggplot(x, aes(box, foldchange.log2, fill=classification)) + geom_violin() + coord_flip())

  # trim the foldchanges for violin plot
  x$FC = x$foldchange.log2
  for(b in x$box) {
    rows = x$box == b
    b_lims = quantile(x$foldchange.log2[rows], c(.01,.99), na.rm=T)
    x$FC[rows & (x$FC < b_lims[1] | x$FC > b_lims[2])] = NA
  }
  print(ggplot(x, aes(box, FC, fill=classification)) +
          geom_violin(na.rm=T) + # draw_quantiles = 0.5   =  deceptive, as we already remove 'outliers' above
          scale_x_discrete(labels=sub(".* ", "", levels(x$box))) +
          labs(y="log2 foldchange between groups in contrast", x="", title = id) +
          coord_flip() +
          theme_bw() + theme(plot.title = element_text(size=8)))
}



#### regular plots; volcano
plain_volcano = function(tib, mtitle="") {
  xlim_symmetric = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(0.001, 0.999), na.rm=T)))
  print(ggplot(tib, aes(foldchange.log2, qvalue_minlog10, colour = color_code, fill = color_code)) +
          geom_point(alpha = 0.5, na.rm=T) +
          scale_color_brewer(palette="Dark2", aesthetics = c("colour", "fill"), guide = guide_legend(title = "", direction = "horizontal", nrow = 1)) +
          facet_wrap( ~ algo_de, ncol = 2, scales = "free_y") +
          xlim(xlim_symmetric) +
          labs(x = "log2 foldchange (far extremes not shown)", y = "-log10 FDR adjusted p-value", title = mtitle) +
          theme_bw() +
          theme(plot.title = element_text(size=8), legend.position = "bottom"))
}
for(id in unique(roc_data$dataset_contrast_norm_algode_id)) { # id=unique(roc_data$dataset_contrast_norm_algode_id)[1]
  x = roc_data %>%
    filter(dataset_contrast_norm_algode_id == id & classification %in% c("background", "foreground")) %>%
    mutate(color_code = paste0(classification, " signif:", signif))
  plain_volcano(x %>% filter(classification == "background"), id)
  plain_volcano(x, id)
}



#### regular plots; compare algo_de within same normalization
for(id in unique(roc_data$dataset_contrast_id)) {
  for(alg in unique(roc_data$algo_norm)) { #id=roc_data$dataset_contrast_id[1]; alg="vwmb"
    x = tib_roc %>% filter(dataset_contrast_id == id & algo_norm == alg) %>% arrange(algo_de)
    for(i in 1:nrow(x)) {
      pROC::plot.roc(x$proc_pauc[[i]], add=(i>1), col=i, main=paste(id, alg, "|=youden \\=signif"), lwd=1, cex.main=0.6)
      points(x$best_specificity[i], x$best_sensitivity[i], col=i, pch="|", cex=1.5)
      points(x$signif_specificity[i], x$signif_sensitivity[i], col=i, pch="\\", cex=1.5)
    }
    legend("bottomright", legend = sprintf("pAUC:%.2f%% TP:%s FP:%s TPR:%.2g FPR:%.2g %s", x$pauc*100, x$signif_tp, x$signif_fp, x$signif_tpr, x$signif_fpr, x$algo_de), col=1:nrow(x), bty='n', cex=.6, lty=1)
  }
}

#### regular plots; compare normalization within same algo_de
for(id in unique(roc_data$dataset_contrast_id)) {
  for(alg in unique(roc_data$algo_de)) { #id=roc_data$dataset_contrast_id[1]; alg="ebayes"
    x = tib_roc %>% filter(dataset_contrast_id == id & algo_de == alg) %>% arrange(desc(pauc))
    for(i in 1:nrow(x)) {
      pROC::plot.roc(x$proc_pauc[[i]], add=(i>1), col=i, main=paste(id, alg), lwd=1, cex.main=0.6)
    }
    legend("bottomright", legend = sprintf("pAUC:%.2f%% TP:%s FP:%s TPR:%.2g FPR:%.2g %s", x$pauc*100, x$signif_tp, x$signif_fp, x$signif_tpr, x$signif_fpr, x$algo_norm), col=1:nrow(x), bty='n', cex=.6, lty=1)
  }
}

dev.off()





pdf("C:/temp/pauc_between_normalizations.pdf", width = 3, height = 3, pointsize = 11, family = "Helvetica")

####
x = tib_roc %>% group_by(dataset_contrast_id, algo_de) %>% mutate(pauc_dist = (pauc - max(pauc)) * 100)
# x = tib_roc %>% group_by(dataset_contrast_id) %>% mutate(pauc_dist = (pauc - max(pauc)) * 100) # obsoleted, compare between algo_de as well
# sort normalization algorithms by overall performance
# TODO: after ggplot release, update to use plot.title.position  @  https://github.com/tidyverse/ggplot2/pull/3494
ref_algo_order = x %>% group_by(algo_norm) %>% summarise(pauc_dist=mean(pauc_dist)) %>% arrange(desc(pauc_dist))
x$algo_norm = factor(x$algo_norm, levels=ref_algo_order$algo_norm)
print(ggplot(x, aes(pauc_dist, algo_norm)) +
        geom_jitter(size = 1.2, height = 0.1, shape = 1, stroke = 0.3) +
        # geom_point(size = 0.5) +
        geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red") +
        labs(caption="compare pAUC between normalizations:\nwithin each DE algorithm in each datasets*contrasts, compare pAUC to the normalization with 'best' pAUC\nenables between-normalization comparison, since performance of DE algorithms differs within and between datasets", x="pAUC distance from max within same dataset", y="") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0, size=6),
              axis.title.x = element_text(size=8, hjust = 1),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8),
              plot.caption = element_text(size=4)) )

dev.off()


# plot separate by algo_de
plotlist_by_algode = plotlist_by_algode_datasetcontrast = list()
for(alg in unique(tib_roc$algo_de)) { # alg="ebayes"
  # subset plot data
  y = x %>% filter(algo_de == alg)
  # update sorting and mean value within current algo_de
  ref_algo_order = y %>% group_by(algo_norm) %>% summarise(pauc_dist=mean(pauc_dist)) %>% arrange(desc(pauc_dist))
  levels(y$algo_norm) = ref_algo_order$algo_norm

  # print(ggplot(y, aes(pauc_dist, algo_norm)) + geom_point() + geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red") +
  #         labs(title=sprintf("consider only %s; compare pAUC between normalizations within the same datasets*contrasts", alg), subtitle = alg, x="pAUC distance from max within same dataset", y="") + theme_bw() + theme(plot.title = element_text(size=8)))
  # print(ggplot(y, aes(pauc_dist, algo_norm, colour=dataset_contrast_id)) + geom_point() + geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red") +
  #         labs(title=sprintf("consider only %s; compare pAUC between normalizations within the same datasets*contrasts", alg), subtitle = alg, x="pAUC distance from max within same dataset", y="") + theme_bw() + theme(plot.title = element_text(size=8), legend.position = "bottom", legend.direction = "vertical"))

  p = ggplot(y, aes(pauc_dist, algo_norm)) +
    # geom_point(size = 0.5) +
    # geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red") +
    labs(caption=sprintf("consider only %s; compare pAUC between normalizations within the same datasets*contrasts", alg), x="pAUC distance from max within same dataset", y="") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0, size=6),
          axis.title.x = element_text(size=8, hjust = 1),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8),
          plot.caption = element_text(size=4),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size=4))

  plotlist_by_algode[[alg]] = p + geom_jitter(size = 1.2, height = 0.1, shape = 1, stroke = 0.3) +
    geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red")
  plotlist_by_algode_datasetcontrast[[alg]] = p + geom_jitter(aes(colour=dataset_contrast_id), size = 1.2, height = 0.1, alpha = 0.6) +
    geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red")
}

pdf("C:/temp/pauc_between_normalizations_by_algo_de.pdf", width = 3, height = 3, pointsize = 11, family = "Helvetica")
for(p in plotlist_by_algode) print(p)
dev.off()
pdf("C:/temp/pauc_between_normalizations_by_algo_de_colorcode_dataset.pdf", width = 5, height = 7, pointsize = 11, family = "Helvetica")
for(p in plotlist_by_algode_datasetcontrast) print(p)
dev.off()








pdf("C:/temp/new_pauc_tpr_fpr__combined_metric.pdf", width = 7, height = 3.5, pointsize = 11, family = "Helvetica")

#### define set of significant proteins at each dataset*contrast*norm; proteins significant in at least 2 DE algorithms (considering ebayes, msqrob, msempire)
MS_prot_signif = roc_data %>%
  # filter significant proteins and algo_de of choice
  filter(signif & algo_de %in% c("ebayes", "msqrob", "msempire")) %>%
  group_by(dataset_contrast_norm_id, protein_id, classification) %>%
  # count how often each protein_id is listed
  count(protein_id, name = "signif_count") %>%
  # select those signif at least 2 times
  filter(signif_count >= 2) %>%
  # for whatever remains with new significance definition, tally the classifications per dataset*contrast*normalization (eg; how many background, how many foreground)
  group_by(dataset_contrast_norm_id) %>%
  count(classification, name = "n_signif")

# total number of true/false positives  +  dataset metadata
MS_classification_by_contrast_norm = roc_data %>%
  distinct(dataset_contrast_norm_id, protein_id, .keep_all = T) %>%
  count(dataset_contrast_norm_id, classification, name = "n_total")

# summarize TP/FP and respective rates
merged_significant = left_join(MS_prot_signif, MS_classification_by_contrast_norm) %>%
  pivot_wider(id_cols = dataset_contrast_norm_id, names_from = classification, values_from = c(n_signif, n_total)) %>%
  group_by(dataset_contrast_norm_id) %>%
  summarise(fp = n_signif_background,
            tp = n_signif_foreground,
            fpr = n_signif_background / n_total_background,
            tpr = n_signif_foreground / n_total_foreground) %>%
  left_join(roc_data %>% distinct(dataset_contrast_norm_id, algo_norm, dataset_name, contrast))


####
for(alg in unique(tib_roc$algo_norm)) { #alg="vwmb"
  # get the pAUC, FPR and TPR for current normalization algorithm as computed by each algo_de
  # then add values from our combined approach
  x = bind_rows(tib_roc %>% filter(algo_norm == alg) %>%
                  mutate(fpr = (1 - signif_specificity),
                         tpr = signif_sensitivity) %>%
                  select(dataset_name, contrast, pauc, tpr, fpr, algo_de),
                #
                merged_significant %>% filter(algo_norm == alg) %>%
                  mutate(pauc = NA, algo_de = "*combined") %>%
                  select(dataset_name, contrast, pauc, tpr, fpr, algo_de))

  ## prepare for plot
  # well formatted grouping variable (dataset*contrast)
  x$id = paste(x$dataset_name, x$contrast)
  x$id = factor(x$id, levels = as.vector(outer(sort(unique(x$dataset_name)), sort(unique(x$contrast)), paste)) )
  x$id = droplevels(x$id)
  # fixed order for algo_de
  x$algo_de = factor(x$algo_de, levels = c("ebayes", "msqrob", "msempire", "msqrobsum", "*combined"))
  # coordinate; position in factor(level) + some minor offset based on algo_de factor(level)
  x$id_coord = as.numeric(x$id) + (as.numeric(x$algo_de) / length(levels(x$algo_de)) - 0.5) / 5  # as.numeric(<factor>) is the same as match(x$id, levels(x$id))
  # convert to percentages
  x$fpr = x$fpr * 100
  x$tpr = x$tpr * 100
  x$pauc = x$pauc * 100

  # for the actual plot, we convert the wide-format data to long-format so we can use facets
  tib_plot = x %>%
    select(id, id_coord, algo_de, fpr, tpr, pauc) %>%
    pivot_longer(cols = c(fpr, tpr, pauc), names_to = "name", values_to = "value")
  tib_plot$name = factor(tib_plot$name, levels = c("pauc", "tpr", "fpr"))

  # plot
  p = ggplot(tib_plot, aes(value, id_coord, colour=algo_de, shape=algo_de)) +
    geom_point(na.rm=T) +
    scale_y_continuous(breaks = seq_along(levels(tib_plot$id)), labels=levels(tib_plot$id)) +
    facet_grid(cols = vars(name), scales = "free_x",
               labeller = labeller(name = array(c("pAUC", "True Positive Rate", "False Positive Rate"), dimnames = list(c("pauc", "tpr", "fpr"))) ) ) +
    theme_bw() +
    scale_color_discrete(guide = guide_legend(title = "Differential Expression Algorithm", direction = "horizontal", title.position = "top", nrow = 2)) +
    scale_shape(guide = guide_legend(title = "Differential Expression Algorithm", direction = "horizontal", title.position = "top", nrow = 2)) +
    labs(title = sprintf("normalization: %s", alg), x = "", y = "") +
    theme(plot.title = element_text(size=8),
          legend.position = "bottom",
          # legend.title = element_blank(),
          # legend.direction = "horizontal",
          legend.title = element_text(size=8),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=6))

  print(p)
}

dev.off()






















########################################################################################################################################################################
#### paper figure
plotlist = list()

#### ROC

# unique(roc_data$dataset_contrast_id)
x = tib_roc %>% filter(dataset_contrast_id == "lfqbench_spectronaut13-nonorm_mindetect3 contrast: A vs B" & algo_de == "ebayes") %>% mutate(algo_norm = sub("&modebetween", "+MB", algo_norm, fixed=T))
# x_clr = RColorBrewer::brewer.pal(n=nrow(x), "Set1") #rainbow(nrow(x), alpha = 0.8)

# hardcoded: color-code each normalization algorithm
x$algo_norm_base = sub("\\+.*", "", x$algo_norm)
# x$algo_norm_base = sub("&.*", "", x$algo_norm)
clr = c("vwmb"="#E41A1C", "msempire"="#377EB8", "loess"="#4DAF4A", "rlr"="#984EA3", "vsn"="#FF7F00")
ggplot_clr = array(clr[match(x$algo_norm_base, names(clr))], dimnames = list(x$algo_norm))[order(match(x$algo_norm_base, names(clr)))]

x = x %>% arrange(match(algo_norm_base, names(clr)))

# hardcoded: normalization algorithms that we want in each plot + hardcoded order
proc_list = x$proc
names(proc_list) = x$algo_norm
proc_list_plot = list(A=proc_list[match(c("vwmb", "msempire", "rlr", "vsn", "loess"), names(proc_list))],
                      C=proc_list[match(c("vwmb", "msempire", "rlr+MB", "vsn+MB", "loess+MB"), names(proc_list))])
for(n in names(proc_list_plot)) { #n="A"
  l = rev(proc_list_plot[[n]]) # reverse so we plot the 'best' methods last
  plotlist[[n]] = pROC::ggroc(l, aes="colour") +
    geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed") +
    scale_colour_manual(values=ggplot_clr, breaks = rev(names(l))) + #
    # labs(x="FPR", y="TPR") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position=c(1, 0),
          legend.justification = c(1, 0),
          legend.spacing.y = unit(0.05, 'lines'),
          legend.key.height = unit(0.05, 'lines'),
          legend.key.width = unit(0.5, 'lines'),
          legend.background = element_blank(),
          # legend.margin = unit(c(2,2,2,2), "pt"),
          legend.text = element_text(size = rel(0.7)),
          legend.title = element_blank())
}


#### volcano

plot_volcano_2 = function(tib) {
  xlim_symmetric = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(0.001, 0.999), na.rm=T)))
  ggplot(tib, aes(foldchange.log2, qvalue_minlog10, colour = color_code, fill = color_code)) +
    geom_hline(yintercept = 2, colour="grey", linetype="dashed", lwd=1) +
    geom_point(alpha = 0.5, na.rm=T) +
    scale_color_brewer(palette="Dark2", aesthetics = c("colour", "fill"), guide = guide_legend(title = "", direction = "horizontal", nrow = 1)) +
    xlim(xlim_symmetric) +
    labs(x = "log2 foldchange", y = "-log10 q-value") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

plotlist$B = plot_volcano_2(roc_data %>%
                              filter(dataset_contrast_norm_algode_id == "lfqbench_spectronaut13-nonorm_mindetect3 contrast: A vs B ebayes @ vsn" & classification == "background") %>%
                              mutate(color_code = paste0(classification, " signif:", signif)))
plotlist$D = plot_volcano_2(roc_data %>%
                              filter(dataset_contrast_norm_algode_id == "lfqbench_spectronaut13-nonorm_mindetect3 contrast: A vs B ebayes @ vsn&modebetween" & classification == "background") %>%
                              mutate(color_code = paste0(classification, " signif:", signif)))

#### pauc scores, copy/paste from pauc_between_normalizations.pdf code
x = tib_roc %>% mutate(algo_norm = sub("&modebetween", "+MB", algo_norm, fixed=T)) %>% group_by(dataset_contrast_id, algo_de) %>% mutate(pauc_dist = (pauc - max(pauc)) * 100)
# sort normalization algorithms by overall performance
ref_algo_order = x %>% group_by(algo_norm) %>% summarise(pauc_dist=mean(pauc_dist)) %>% arrange(desc(pauc_dist))
x$algo_norm = factor(x$algo_norm, levels=ref_algo_order$algo_norm)

plotlist$E = ggplot(x, aes(pauc_dist, algo_norm)) +
  geom_jitter(size = 1.2, height = 0.1, shape = 1, stroke = 0.3, na.rm = T) +
  # geom_point(size = 0.5) +
  geom_point(data=ref_algo_order, aes(pauc_dist, algo_norm), shape = "|", size=5, colour="red") +
  labs(x="pAUC distance from optimum", y="") +
  theme_bw() +
  theme(axis.title.x.bottom = element_text(hjust = 1))



#### combine
l = plotlist[LETTERS[1:4]]
p_comb_abcd = ggpubr::ggarrange(plotlist = l, ncol = 2, nrow=2, labels=names(l), font.label = list(size = 11, color = "black", face = "bold", family = NULL))
p_comb_abcde = ggpubr::ggarrange(p_comb_abcd, plotlist$E, ncol=2, nrow=1, widths = c(2,1), labels=c("", "E"), font.label = list(size = 11, color = "black", face = "bold", family = NULL))

graphics.off(); Sys.sleep(3) # bugfix to force previous graphics device to close
pdf("C:/temp/fig2.pdf", width = 6, height = 4, pointsize = 11, family = "Helvetica") # A4 full size = 8.27 x 11.69 inches
print(p_comb_abcde)
dev.off()

########################################################################################################################################################################










# print(ggplot(x, aes(pauc * 100, id, colour=algo_de, shape=algo_de)) +
#         geom_point(aes(y=id_coord)) +
#         geom_point(alpha = 0) +
#         scale_y_discrete(name = "", labels = levels(x$id)) +
#         labs(title=paste("ROC pAUC statistic per dataset. normalization:", alg), x="partial area under the curve (in %), at 5% FPR in ROC", y="") +
#         theme_bw() + theme(plot.title = element_text(size=8), legend.position = "bottom", axis.text.x = element_text(size=8), axis.text.y = element_text(size=6)))
#
#
# xx = x %>% mutate(fpr = (1-signif_specificity), tpr = signif_sensitivity) %>% select(id, tpr, fpr, algo_de)
# # test; add 'merged DE'
# xx_mergedsignif = merged_significant %>% filter(algo_norm == alg)
# xx_mergedsignif$id = paste(xx_mergedsignif$dataset_name, xx_mergedsignif$contrast)
# xx_mergedsignif$algo_de = "*combined"
# xx = bind_rows(xx, xx_mergedsignif %>% select(id, tpr, fpr, algo_de))
# xx$algo_de = factor(xx$algo_de, levels = levels(x$algo_de))
# xx$id = factor(xx$id, levels = levels(x$id))
# xx$id_coord = as.numeric(xx$id) + (as.numeric(xx$algo_de) / length(levels(xx$algo_de)) - 0.5) / 5
#
# ## plot FPR at significance cutoff
# print( ggplot(xx, aes(fpr * 100, id, colour=algo_de, shape=algo_de)) +
#          geom_point(aes(y=id_coord)) +
#          geom_point(alpha = 0) +
#          scale_y_discrete(name = "", labels = levels(xx$id)) +
#          labs(title=paste("ROC per dataset; how conservative is each algo_de? normalization:", alg), x="False Positive Rate (in %) at qvalue<=0.01", y="") +
#          theme_bw() + theme(plot.title = element_text(size=8), legend.position = "bottom", axis.text.x = element_text(size=8), axis.text.y = element_text(size=6)) )
# # for code QC, minimal plot: ggplot(xx, aes(fpr * 100, id, colour=algo_de, shape=algo_de)) + geom_point()
#
# # copy/paste of above, but for TPR
# print( ggplot(xx, aes(tpr * 100, id, colour=algo_de, shape=algo_de)) +
#          geom_point(aes(y=id_coord)) +
#          geom_point(alpha = 0) +
#          scale_y_discrete(name = "", labels = levels(xx$id)) +
#          labs(title=paste("ROC per dataset; how conservative is each algo_de? normalization:", alg), x="True Positive Rate (in %) at qvalue<=0.01", y="") +
#          theme_bw() + theme(plot.title = element_text(size=8), legend.position = "bottom", axis.text.x = element_text(size=8), axis.text.y = element_text(size=6)))
#
# # print(ggplot(x %>% arrange(dataset_contrast_id, algo_de), aes(best_specificity, dataset_contrast_id, colour=algo_de, shape=algo_de)) + geom_point() + labs(title=alg) )



#################### test code for benchmark plots
#
# ## plot;
# # load("C:/temp/lfqbench-modified_spectronaut13__detect=3_quant=3_minpep=1_topn=10.RData")
# load("C:/temp/lfqbench_spectronaut13__detect=3_quant=3_minpep=1_topn=10.RData")
# library(pROC)
# D = dataset[[1]]
#
# for(contr in unique(D$stats_de$contrast)) {
#   lbl = paste(contr, "norm:", paste(D$norm_algorithm, collapse="&"))
#   stats_contr = D$stats_de %>% filter(contrast == contr)
#
#   print(plot_foldchanges(stats_contr %>% mutate(color_code = classification), mtitle=lbl))
#
#   tib_roc = stats_contr %>%
#     add_column(predictor = stats_contr$classification) %>%
#     filter(stats_contr$classification %in% c("background", "foreground"))
#
#   plot_roc(tib = tib_roc, universe = "all", mtitle=lbl, plot_coords = TRUE)
# }

