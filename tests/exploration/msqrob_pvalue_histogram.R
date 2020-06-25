# cl <<- initialize_multiprocessing()
#
# ### adapted from https://github.com/statOmics/MSqRobSum/blob/master/vignettes/msqrobsum.Rmd
#
# library(limma)
# library(MSnbase)
#
# # read example dataset @ MSqRob vignette
# data_path = "C:/Users/Frank/Downloads/peptides.txt.gz"
# data_path = gzfile(data_path)
# exprs_col = grepEcols(data_path, 'Intensity ',split = '\t')
# set = readMSnSet2(data_path ,ecol = exprs_col,fnames = 'Sequence', sep = '\t',stringsAsFactors = FALSE)
#
# # blindly remove any protein ID that indicates reverse/contaminant
# set <- set[!grepl("__", fData(set)$Proteins)]
#
# # add sample metadata
# sampleNames(set) = str_replace(sampleNames(set),'Intensity.','')
# pd = data.frame(condition = as.factor(str_extract(sampleNames(set),'^.')))
# rownames(pd) = sampleNames(set)
# pData(set) = pd
#
# # remove NA
# exprs(set)[0 == (exprs(set))] <- NA
# # plotNA(set)
#
# # normalize
# # plotMDS(exprs(log(set, base = 2)), top = Inf,col = as.integer(pData(set)$condition))
# set = MSnbase::normalize(set, 'vsn')
# # plotMDS(exprs(set), top = Inf,col = as.integer(pData(set)$condition))
#
# # msqrob
# formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
# msqrob_result <- msqrobsum(data = set, formulas, contrasts = 'condition', mode = 'msqrob', group_vars = c('Proteins'))
# msqrob_result
#
# contrasts = msqrob_result %>%
#   select(Proteins, contrasts) %>%
#   unnest(cols = contrasts)
# filter(contrasts,qvalue <.05) %>% group_by(contrast) %>% tally()
#
#
#
# ### finally, plot a p-value histogram
# hist(contrasts %>% filter(contrast == "conditionb-conditiona") %>% pull(pvalue), breaks = 20, main = "https://github.com/statOmics/MSqRobSum/blob/master/vignettes/msqrobsum.Rmd", cex.main=0.8)
# hist(contrasts %>% filter(contrast == "conditionb-conditiona") %>% pull(qvalue), breaks = 20)
#
