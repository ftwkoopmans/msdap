testthat::context("normalization reproducibility")
msdap::enable_log(FALSE)
datasets = prepare_test_datasets(5000)
datasets_as_matrix = datasets$as_matrix

for(lbl in names(datasets_as_matrix)) {
  m_group_id = datasets_as_matrix[[lbl]]$groupid
  m_log2 = datasets_as_matrix[[lbl]]$mat
  m_log2[!is.finite(m_log2)] = NA
  # cat("\ndataset:", lbl, "\n")

  for(alg_norm in msdap::normalization_algorithms()) {
    # cat(alg_norm, "...", sep ="")
    # try N times + ensure we don't rely on some random seed state by generating a random number in between
    m_log2_norm_1 = msdap::normalize_matrix(m_log2, algorithm = alg_norm, group_by_cols = m_group_id)
    testthat::test_that(paste("dataset:", lbl, "norm:", alg_norm), {
      for(i in 1:3) {
        # cat(" ", i)
        rng = rnorm(1)
        m_log2_norm_2 = msdap::normalize_matrix(m_log2, algorithm = alg_norm, group_by_cols = m_group_id)
        testthat::expect_equal(m_log2_norm_1, m_log2_norm_2)
        # i_is_eq = all.equal(m_log2_norm_1, m_log2_norm_2)
        # stopifnot(i_is_eq)
        # plot(m_log2_norm_1, m_log2_norm_2, main = paste(alg_norm, i))
      }
    })
    # cat("\n")
  }
}
