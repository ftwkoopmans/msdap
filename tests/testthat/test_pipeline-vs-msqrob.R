testthat::context("assert pipeline output equal to example from MSqRobSum R package")
msdap::enable_log(FALSE)
cl <<- msdap::initialize_multiprocessing()


######################################################## generate results ########################################################

library(limma)
library(MSnbase)
library(nloptr)

### following example code from MSqRobSum vignette; read example dataset included with MSqRobSum package
# https://github.com/statOmics/MSqRobSum/blob/master/vignettes/msqrobsum.Rmd
data_path = system.file('extdata','peptides.txt.gz', package = 'msqrobsum')
data_path = gzfile(data_path)
exprs_col = MSnbase::grepEcols(data_path, 'Intensity ',split = '\t')
set = MSnbase::readMSnSet2(data_path ,ecol = exprs_col,fnames = 'Sequence', sep = '\t',stringsAsFactors = FALSE)
# remove redundant words in sample names  &  add sample metadata
MSnbase::sampleNames(set) = stringr::str_replace(MSnbase::sampleNames(set),'Intensity.','')
pd = data.frame(condition = as.factor(stringr::str_extract(MSnbase::sampleNames(set),'^.')))
rownames(pd) = MSnbase::sampleNames(set)
MSnbase::pData(set) = pd
### end code snippet copied from documentation

### custom filtering step, only needed for this input file because it is some custom edit of the original MaxQuant peptides.txt file
### it appears this dataset still contains contaminants but not the respective column
set <- set[!grepl("__", MSnbase::fData(set)$Proteins) & MSnbase::fData(set)$Reverse == "", ]

# remove NA + log transform
## my code, same result
x = MSnbase::exprs(set)
x[!is.finite(x) | x == 0] = NA
Biobase::exprs(set) <- log2(x)

### custom data filter; only conditions A and B, only peptides that have 2+ values in each 'group of samples'
rows = rowSums(!is.na(MSnbase::exprs(set)[,pd$condition=="a"])) >= 2 &
  rowSums(!is.na(MSnbase::exprs(set)[,pd$condition=="b"])) >= 2
set_filtered = set[rows, pd$condition %in% c("a", "b")]

# optionally normalize
# set_filtered = MSnbase::normalize(set_filtered, 'vsn')

# msqrob
formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
msqrob_output = msqrobsum(data = set_filtered, formulas, contrasts = 'condition', mode = 'msqrob', group_vars = c('Proteins'))
msqrob_result = msqrob_output %>% select(Proteins, contrasts) %>% unnest(cols = contrasts)

# # print number of significant proteins
# print( msqrob_result %>% filter(qvalue <= 0.01) %>% group_by(contrast) %>% tally() %>% arrange(contrast, n) )



############## analogous for msdap


# parse input data
dataset = msdap::import_maxquant_peptides(file_peptides = system.file('extdata','peptides.txt.gz', package = 'msqrobsum'), remove_shared_peptides = F, pep2prot = NA)
### custom filtering step, only needed for this input file because it is some custom edit of the original MaxQuant peptides.txt file
### it appears this dataset still contains contaminants but not the respective column
dataset$peptides = dataset$peptides %>% dplyr::filter(!grepl("__", protein_id)) # note; reverse hits were already removed while importing peptides.txt

# extract groups from filenames by specifying a regex for characters to remove from the sample name (whatever is not removed = groupn name)
dataset = msdap::sample_metadata_custom(dataset, sample_property_regex = list(shortname = "", group = "\\d")) #, sample_exclude_regex = "^(c|d)")
dataset = msdap::setup_contrasts(dataset, contrast_list = list(c("a", "b"))) # setup a contrast, just compare groups A and B
# dataset$samples

# filter N detect in both groups
dataset = msdap::filter_dataset(dataset,
                                filter_min_detect = 2,
                                norm_algorithm = "",
                                by_group = F,
                                all_group = F,
                                by_contrast = T)

# DE analysis
dataset = msdap::dea(dataset, qval_signif = 0.01, algo_de = "msqrob") # output_dir_for_eset = "C:/temp/")
msdap_result = dataset$de_proteins
# print( msdap_result %>% filter(signif) %>% group_by(contrast, algo_de) %>% tally() %>% arrange(contrast, n) )



############## peptide abundances are the same throughout code bases


### there is an incredibly minor imprecision between intensities generated in our pipeline versus above MSqRob example code
### this leads to minor differences downstream, a technicality that cannot be helped (not found a workaround yet)
### note; if you replace the intensities used for msdap:dea() upstream the p-values are 100% the same

# take ExpressionSet from MSqRobSum code above, and convert the log2 peptide intensity matrix to long format tibble
tmp = dplyr::as_tibble(MSnbase::exprs(set_filtered)) %>%
  tibble::add_column(peptide_id=rownames(MSnbase::exprs(set_filtered))) %>%
  tidyr::pivot_longer(cols = c(-peptide_id), names_to = "sample_id", values_to = "intensity") %>%
  dplyr::rename(intensity_eset = intensity)
# join with peptides tibble and compare intensity values
peptides = dplyr::left_join(dataset$peptides, tmp) #set this if you want to use the peptide tibble in msdap downstream; attr(peptides, 'acquisition_mode') <- 'dda'

# assert similarity
testthat::test_that("peptide abundance are the same between msdap and msqrobsum code ?",
                    testthat::expect_equal(peptides$`intensity_contrast: a vs b`, peptides$intensity_eset, tolerance = sqrt(.Machine$double.eps)))
print("numerical imprecision in peptide intensity differences")
print(all.equal(peptides$`intensity_contrast: a vs b`, peptides$intensity_eset, tolerance = 0)) # see the difference if we require zero tolerance

# print("peptide abundance are the same between msdap and msqrobsum code ?")
# print(all.equal(peptides$`intensity_contrast: a vs b`, peptides$intensity_eset, tolerance = sqrt(.Machine$double.eps))) # default tolerance
# print("but is it EXACTLY the same?")
# print(all.equal(peptides$`intensity_contrast: a vs b`, peptides$intensity_eset, tolerance = 0)) # see the difference if we require zero tolerance




############## p-values are the same

# data tables are aligned
testthat::test_that("same dimensions", testthat::expect_equal(nrow(msqrob_result), nrow(msdap_result)))
testthat::test_that("protein id order", testthat::expect_equal(msqrob_result$Proteins, msdap_result$protein_id))
testthat::test_that("log2 fold-change", testthat::expect_equal(msqrob_result$logFC, msdap_result$foldchange.log2, tolerance = 10^-3))
testthat::test_that("pvalue", testthat::expect_equal(msqrob_result$pvalue, msdap_result$pvalue, tolerance = 10^-3))
print("exact p-value equality due to numerical imprecision:")
print(all.equal(msqrob_result$pvalue, msdap_result$pvalue))

# plot(log10(msqrob_result$pvalue), log10(msdap_result$pvalue))
# plot(msqrob_result$logFC, msdap_result$foldchange.log2)


# finally, shut down multithreading clusters
parallel::stopCluster(cl)







# sample_of_interest = samples %>% filter(group %in% c("a", "b")) %>% pull(sample_id)
# msdap_eset = msdap::tibble_as_eset(peptides %>% filter(is.finite(`intensity_contrast: a vs b`) & sample_id %in% sample_of_interest), proteins, samples)
#
# # test: msdap_eset -->> insert directly into msdap::de_msempire() -->> assert same as above
# pData(msdap_eset)$condition = pData(msdap_eset)$group
# formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
# msdap_eset_msqrob_output = msqrobsum(data = MSnbase::as.MSnSet.ExpressionSet(msdap_eset), formulas, contrasts = 'condition', mode = 'msqrob', group_vars = c('protein_id'))
# msdap_eset_msqrob_result = msdap_eset_msqrob_output %>% select(protein_id, contrasts) %>% unnest(cols = contrasts)
# # head(pData(msdap_eset)); head(fData(msdap_eset))
# all.equal(msqrob_result$pvalue, msdap_eset_msqrob_result$pvalue)
# all.equal(msdap_eset_msqrob_result$pvalue, msdap_result$pvalue)
#
#
#
# #################
#
# # if this works, filter function introduces rng
# msdap_eset = msdap::tibble_as_eset(peptides %>% filter(is.finite(intensity) & sample_id %in% colnames(exprs(set_filtered)) & peptide_id %in% rownames(exprs(set_filtered))), proteins, samples)
#
# # test: msdap_eset -->> insert directly into msdap::de_msempire() -->> assert same as above
# pData(msdap_eset)$condition = pData(msdap_eset)$group
# formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
# msdap_eset_msqrob_output = msqrobsum(data = MSnbase::as.MSnSet.ExpressionSet(msdap_eset), formulas, contrasts = 'condition', mode = 'msqrob', group_vars = c('protein_id'))
# msdap_eset_msqrob_result = msdap_eset_msqrob_output %>% select(protein_id, contrasts) %>% unnest(cols = contrasts)
# all.equal(msqrob_result$pvalue, msdap_eset_msqrob_result$pvalue)
#
#
# #################
#
# formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
# msqrob_output = msqrobsum(data = set_filtered, formulas, contrasts = 'condition', mode = 'msqrob', group_vars = c('Proteins'))
# head(pData(set_filtered)); head(fData(set_filtered))
#
#
# stopifnot(ncol(set_filtered) == ncol(msdap_eset))
# stopifnot(nrow(set_filtered) == nrow(msdap_eset))
# stopifnot(colnames(set_filtered) == colnames(msdap_eset))
# stopifnot(rownames(set_filtered) == rownames(msdap_eset))
#
# # confirm that the exact same peptide intensities are preserved through either codebase
# stopifnot(sum(!is.na(exprs(set_filtered))) == sum(!is.na(exprs(msdap_eset))))
# all.equal(c(exprs(set_filtered)), c(exprs(msdap_eset)))
# # plot(c(exprs(set_filtered)), c(exprs(msdap_eset)))
# # stopifnot(na.omit(c(exprs(set_filtered))) == na.omit(c(exprs(msdap_eset))))
#
#
# # 2) statistical results are the same
#
# stopifnot(nrow(msqrob_result) == nrow(msdap_result))
# stopifnot(msqrob_result$Proteins == msdap_result$protein_id)
#
# all.equal(msqrob_result$logFC, msdap_result$foldchange.log2)
# all.equal(msqrob_result$pvalue, msdap_result$pvalue)
# plot(log10(msqrob_result$pvalue), log10(msdap_result$pvalue), main = "msqrob vs msdap: log10(pvalue)")
#
#
# #
# #
# # ### now convert the same data to a format compatible with our pipeline
# # l = import_expressionset(msempire_eset_filtered)
# #
# # proteins = l$proteins
# # peptides = l$peptides
# # samples = sample_contrasts(l$samples, contrast_list = list(c("0", "1")))
# #
# # # don't strictly need this, but just include all standard operations to ensure results are the same
# # # all filters; global, local, by_group (for CoV plots). adds additional columns to peptide table with filtered and normalized peptide intensities
# # # peptides = filter_dataset(peptides, proteins, samples,
# # #                           filter_min_detect = 2,
# # #                           normalization = "",
# # #                           filter_by_group = F, # if we are making a QC report downstream, prep this data. otherwise, don't bother
# # #                           filter_all_groups = F,
# # #                           filter_by_contrast = T)
# # # DEA
# # msdap_result = dea(peptides, proteins, samples, qval_signif = 0.01, do_ebayes = F, do_msempire = T, do_msqrob = F, do_msqrobsum = F) # output_dir_for_eset = "C:/temp/"
# # msdap_result2 = dea(peptides, proteins, samples, qval_signif = 0.01, do_ebayes = F, do_msempire = T, do_msqrob = F, do_msqrobsum = F) # output_dir_for_eset = "C:/temp/"
# # stopifnot(msdap_result$pvalue == msdap_result2$pvalue)
# #
# # # print( msdap_result %>% filter(signif) %>% group_by(contrast, algo_de) %>% tally() %>% arrange(contrast, n) ) # print number of significant proteins
# #
# #
# #
# #
# # ######################################################## assert equality ########################################################
# #
# # # TODO: do msqrob unit-test first to ensure msdap parts are fine
# #
# # stopifnot(nrow(msempire_result) == nrow(msdap_result))
# # stopifnot(msempire_result$prot.id == msdap_result$protein_id)
# # stopifnot(msempire_result$log2FC == msdap_result$foldchange.log2)
# # stopifnot(msempire_result$p.val == msdap_result$pvalue)
# #
# # all.equal(msempire_result$p.val, msdap_result$pvalue)
# # all.equal(msempire_result$log2FC, msdap_result$foldchange.log2)
# #
# # testthat::test_that("same dimensions", testthat::expect_equal(nrow(msempire_result), nrow(msdap_result)))
# # testthat::test_that("protein id order", testthat::expect_equal(msempire_result$prot.id, msdap_result$protein_id))
# # testthat::test_that("log2 fold-change", testthat::expect_equal(msempire_result$log2FC, msdap_result$foldchange.log2))
# # testthat::test_that("pvalue", testthat::expect_equal(msempire_result$p.val, msdap_result$pvalue))
# # testthat::test_that("qvalue", testthat::expect_equal(msempire_result$p.adj, msdap_result$qvalue))
# #
