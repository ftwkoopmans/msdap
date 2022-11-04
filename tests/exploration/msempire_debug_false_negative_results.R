# library(msEmpiRe)
# library(MSnbase)
#
# # load ExpressionSet we prepared earlier; LFQbench dataset processed by Spectronaut. we selected only peptides with good q-value and applied vsn normalization
# load("C:/temp/ExpressionSet_contrast A vs B.RData")
#
# # we stored our intensity data on log2 scale, transform back to regular intentities for ms-empire
# x = exprs(eset)
# x = 2^x
# x[!is.finite(x)] = 0
# exprs(eset) = x
# # this eset was filtered previously, there are no missing values
# table(rowSums(x>0))
#
# # we stored our protein identifiers as protein_id, so add the prot.id column expected by ms-empire
# fData(eset)["prot.id"] = fData(eset)[, "protein_id"]
#
# # ms-empire normalize & DE
# eset_norm = msEmpiRe::normalize(eset) #, out.dir = "c:/temp/")
# mse_result = msEmpiRe::de.ana(eset_norm) #, out.dir = "c:/temp/")
#
# # many proteins have a pvalue of zero
# table(mse_result$p.val == 0)
# head(mse_result[mse_result$p.val == 0,])
#
# # volcano shows a strong binning of pvalues, how come? and does it matter?
# # (pvalues that are zero are not shown due to -log10, doesn't matter for this example)
# plot(mse_result$log2FC, -log10(mse_result$p.adj))
#
#
#
# colnames(mse_result)
#
#
# # rollup to protein level, for eBayes we use the traditional sum approach
# prot_int = 2^exprs(eset_norm)
# prot_int[!is.finite(prot_int)] = 0
# eset_tmp = eset_norm
# exprs(eset_tmp) = prot_int
#
# protset = suppressWarnings(suppressMessages(combineFeatures(MSnbase::as.MSnSet.ExpressionSet(eset_tmp), fun = "sum", groupBy = fData(eset)$protein_id))) # fun = ifelse(protein_rollup_robust, "robust", "sum")
# x = log2(exprs(protset))
# x[!is.finite(x)] = NA
# group_by_cols = pData(protset)$condition
# # group_by_cols = match(group_by_cols, unique(group_by_cols)) - 1
#
# # ref implementation: contr_design = model.matrix(~ rep(0:1, c(length(contr_cols_a), length(contr_cols_b))))
# fit = suppressWarnings(eBayes(lmFit(x, model.matrix(~group_by_cols))))
# # !! sort.by="none" keeps the output table aligned with input matrix
# result = suppressMessages(topTable(fit, number = nrow(x), adjust.method = "fdr", sort.by = "none", confint = TRUE))
#
#
# # ###
# #
# #
# #
# #
# # # this protein should be detected as differentially abundant ?
# # # from the data, can only explain this by low abundance + relatively small fold-change. eg; the lower abundance error bin @ msempire strongly penalizes these
# # protein_id = "P06787" # example 1
# # protein_id = "Q12285" # example 2
# #
# #
# # # results from MS-EmpiRe
# # print(result[result$prot.id == protein_id, ])
# # # respective input data from eset
# # y = exprs(eset_norm)[fData(eset_norm)$protein_id == protein_id, ]
# # print(y)
# # # plot, visual inspection
# # plot(seq_along(y), log2(y), las=1)
# #
# #
# # ### find some more examples
# #
# # # variation
# # coefficient_of_variation = function(nonlog_array) sqrt(expm1(sd(log(nonlog_array), na.rm = T)^2))
# # x = exprs(eset_norm)
# # x[x==0] = NA
# # df = data.frame(protein_id = fData(eset_norm)$protein_id,
# #                 cv_group_worst = pmax(apply(x[,1:4], 1, coefficient_of_variation),
# #                       apply(x[,5:8], 1, coefficient_of_variation)),
# #                 fc_est = rowMeans(x[,1:4], na.rm=T) / rowMeans(x[,5:8], na.rm=T) )
# #
# # # protein with single peptide AND abs(fold-change) > 1.5
# # protein_id_signif = result$prot.id[is.finite(result$prot.p.adj) & is.finite(result$prot.p.adj) <= 0.01]
# # prot_count = table(df$protein_id)
# # rows = df$protein_id %in% names(prot_count)[prot_count==1] &
# #   ! df$protein_id %in% protein_id_signif &
# #   abs(log2(df$fc_est)) > log2(1.25)
# #
# # # sort by lowest variation first, apply selection, plot topN
# # df_selection = df[rows, ]
# # df_selection = df_selection[order(df_selection$cv_group_worst, decreasing = F), ]
# #
# # for(pid in head(df_selection$protein_id, 10)) {
# #   y = exprs(eset_norm)[fData(eset_norm)$protein_id == pid, ]
# #   plot(seq_along(y), log2(y), las=1, main=paste("protein_id =", pid))
# # }
# #
