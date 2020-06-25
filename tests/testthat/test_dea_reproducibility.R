testthat::context("DEA reproducibility")
msdap::enable_log(FALSE)
datasets = prepare_test_datasets(1000)
cl <<- msdap::initialize_multiprocessing()

for(lbl in names(datasets$as_tibble)) {
  eset_peptides = datasets$as_eset[[lbl]]

  test_func = list(ebayes = function() msdap::de_interface_ebayes(eset_proteins = msdap::eset_from_peptides_to_proteins(eset_peptides)),
                   msqrob = function() msdap::de_msqrobsum_msqrob(eset_peptides, use_peptide_model = T, input_intensities_are_log2 = T),
                   msempire = function() msdap::de_msempire(eset_peptides, input_intensities_are_log2 = T))
  # cat("\ndataset:", lbl, "\n")
  # head(fData(eset_peptides))
  # head(pData(eset_peptides))
  # hist(exprs(eset_peptides))

  for(algo_de in names(test_func)) { # algo_de = names(test_func)[1]
    res1 = do.call(test_func[[algo_de]], args = list())
    # try N times + ensure we don't rely on some random seed state by generating a random number in between
    testthat::test_that(paste("dataset:", lbl, "DEA:", algo_de), {
      for(i in 1:3) {
        rng = rnorm(1)
        res2 = do.call(test_func[[algo_de]], args = list())
        testthat::expect_equal(res1$protein_id, res2$protein_id)
        testthat::expect_equal(res1$foldchange.log2, res2$foldchange.log2)
        testthat::expect_equal(res1$qvalue, res2$qvalue)
      }
    })
  }
}


# finally, shut down multithreading clusters
parallel::stopCluster(cl)
