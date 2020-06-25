# # Frank Koopmans (frank.koopmans@vu.nl)
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
# # question: many proteins have a pvalue of zero, I assume these are the most significant hits, but why zero ?
# table(mse_result$p.val == 0)
# head(mse_result[mse_result$p.val == 0,])
#
# # question: volcano shows a strong binning of pvalues, how come?
# # (pvalues that are zero are not shown in this plot due to log10 transform, doesn't matter for this example)
# plot(mse_result$log2FC, -log10(mse_result$p.adj))
#
# # same thing, now color-coding by protein species (simple proof of concept, doesn't consider mixed proteinsgroups etc.)
# is_ecoli = grepl("_ECOL", mse_result$prot.id)
# is_yeast = grepl("_YEAS", mse_result$prot.id)
# plot(mse_result$log2FC, -log10(mse_result$p.adj), col = ifelse(is_ecoli, "green", ifelse(is_yeast, "red", "black")))
