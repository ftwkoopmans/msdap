testthat::context("assert pipeline output equal to example from MS-EmpiRe R package")
msdap::enable_log(FALSE)


######################################################## generate results ########################################################

### example dataset and application of MS-EmpiRe following the documentation at: https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R
f <- system.file("extdata", "c1_c3.data", package = "msEmpiRe")
p <- system.file("extdata", "c1_c3.pdata", package = "msEmpiRe")
# loading data from installed data sets. The dataset is the yeast spike-in benchmarking data of O'Connell et al., as discussed in the paper
suppressMessages(data <- msEmpiRe::read.standard(f, p,
                                                 prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                                 signal_pattern="c.*rep.*"))
# extract the first two conditions
conditions <- msEmpiRe::extract_conditions(data)
conditions <- conditions[, c(1,2)]
# removing peptides that are detected in less than 2 samples per condition
tmp = capture.output(msempire_eset_filtered <- msEmpiRe::filter_detection_rate(data, condition=conditions))
### end code snippet copied from documentation



# print("TEST: log conversion rounding errors")
# msempire_eset_filtered_logtest = msempire_eset_filtered
# x = log2(exprs(msempire_eset_filtered))
# x[!is.finite(x)] = NA
# x = 2^x
# x[!is.finite(x)] = 0
# exprs(msempire_eset_filtered_logtest) = x
# set.seed(123)
# capture.output(msempire_result_logtest <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(msempire_eset_filtered_logtest))))
# print("TEST: log conversion rounding errors")



## !! set same seed as in our pipeline before applying MS-EmpiRe
set.seed(123)
capture.output(msempire_result <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(msempire_eset_filtered))))
set.seed(123)
capture.output(msempire_result2 <- suppressWarnings(suppressMessages(msEmpiRe::de.ana(msempire_eset_filtered))))
testthat::test_that("with set.seed(123), msEmpiRe::de.ana() is reproducible", testthat::expect_equal(msempire_result$p.val,  msempire_result2$p.val))



### now convert the same data to a format compatible with our pipeline
dataset = import_expressionset(msempire_eset_filtered)
dataset = setup_contrasts(dataset, contrast_list = list(c("0", "1")))
# DEA
dataset = dea(dataset, qval_signif = 0.01, algo_de = "msempire") # output_dir_for_eset = "C:/temp/"
msdap_result = dataset$de_proteins
# print( msdap_result %>% filter(signif) %>% group_by(contrast, algo_de) %>% tally() %>% arrange(contrast, n) ) # print number of significant proteins


################################
################################
################################
# take ExpressionSet from MSqRobSum code above, and convert the log2 peptide intensity matrix to long format tibble
tmp = dplyr::as_tibble(Biobase::exprs(msempire_eset_filtered)) %>%
  tibble::add_column(peptide_id=rownames(exprs(msempire_eset_filtered))) %>%
  tidyr::pivot_longer(cols = c(-peptide_id), names_to = "sample_id", values_to = "intensity") %>%
  dplyr::rename(intensity_eset = intensity)
tmp$intensity_eset = log2(tmp$intensity_eset)

# join with peptides tibble and compare intensity values
peptides = dataset$peptides %>% dplyr::left_join(tmp)

# assert similarity
testthat::test_that("peptide abundance are the same between msdap and msempire code ?",
                    testthat::expect_equal(peptides$intensity, peptides$intensity_eset, tolerance = 0))
################################
################################
################################



######################################################## assert equality ########################################################

# msempire operates on peptide intensities as integer values, while msdap stores data after log transformation.
# this leads to numerical imprecision, since we convert int->float->int

plot(log10(msempire_result$p.val), log10(msdap_result$pvalue), main = "msempire vs msdap: log10(pvalue)")

testthat::test_that("same dimensions", testthat::expect_equal(nrow(msempire_result), nrow(msdap_result)))
testthat::test_that("protein id order", testthat::expect_equal(msempire_result$prot.id, msdap_result$protein_id))
testthat::test_that("log2 fold-change", testthat::expect_equal(msempire_result$log2FC, msdap_result$foldchange.log2, tolerance = 10^-3))
testthat::test_that("pvalue", testthat::expect_equal(msempire_result$p.val, msdap_result$pvalue, tolerance = 10^-2))

print("exact p-value equality due to numerical imprecision:")
print(all.equal(msempire_result$p.val, msdap_result$pvalue))
# all.equal(msempire_result_logtest$p.val, msdap_result$pvalue)

