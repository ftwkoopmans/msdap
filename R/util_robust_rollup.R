### some test code to compare normal rollup versus our custom rollup
# load(sprintf("C:/DATA/msdap_manuscript_benchmark_datasets/rdata_intermediate_results/%s__pep%s_top%s_detect%s.RData", "OConnel2018_maxquant_mbr", 1, 0, 0))
# tib_dea_ref = dataset$de_proteins %>% filter(is.finite(qvalue) & algo_de == "ebayes") %>% left_join(dataset$proteins %>% mutate(classification = regex_classification(fasta_headers, dataset$protein_classification_regex)) %>% select(protein_id, classification) ) %>% arrange(qvalue)
#
# dataset = dea(dataset, algo_de = "ebayes", rollup_strategy = "ipl")
# tib_dea = dataset$de_proteins %>% filter(is.finite(qvalue) & algo_de == "ebayes") %>% left_join(dataset$proteins %>% mutate(classification = regex_classification(fasta_headers, dataset$protein_classification_regex)) %>% select(protein_id, classification) ) %>% arrange(qvalue)
#
# tib_dea_ref %>% filter(qvalue<=0.01) %>% count(contrast, classification) %>% arrange(contrast, classification)
# tib_dea %>% filter(qvalue<=0.01) %>% count(qvalue<=0.01, contrast, classification) %>% arrange(contrast, classification)


#' @importFrom matrixStats rowMads rowMeans2
#' @importFrom data.table data.table
impute_peptide_local = function(x, protein_ids, groups, par_remove_outliers = TRUE, par_remove_sparse_data = TRUE) {
  list_group_indices = sapply(unique(groups), function(x) which(groups==x), simplify = F)
  stopifnot(lengths(list_group_indices) >= 2) # don't allow groups that only have 1 sample
  # fast aggregation of indices (eg; a table of protein_id -->> list of indices in peptide matrix)
  df_protein_indices = as.data.frame(data.table::data.table(protid=protein_ids, index=seq_along(protein_ids))[, list(index=list(index)), by=protid])

  # expected variation per peptide
  grp_mad = matrix(0, nrow = nrow(x), ncol = length(list_group_indices))
  for(g in seq_along(list_group_indices)) {
    grp_mad[,g] = matrixStats::rowMads(x[,list_group_indices[[g]]], na.rm = T)
  }
  x_mad = matrixStats::rowMeans2(grp_mad, na.rm = T)
  x_mad_mean = mean(x_mad, na.rm = T)

  # init output matrix
  mat_result = matrix(0.0, nrow=nrow(df_protein_indices), ncol=ncol(x), dimnames = list(df_protein_indices$protid, colnames(x)))
  # TODO: use parallel package
  for(protindex in 1:nrow(mat_result)) { #protindex=1
    rows = df_protein_indices$index[[protindex]]
    mat_result[protindex,] = impute_peptide_local_proteinmatrix(x[rows,,drop=F], list_group_indices=list_group_indices, mad_prior = x_mad_mean, par_remove_outliers=par_remove_outliers, par_remove_sparse_data=par_remove_sparse_data)
  }
  return(mat_result)
}



## basic imputation model, assuming pre-filtering such that peptides with very few data points have been removed prior
## value = group mean (or median)  +  sample-specific shift from group-average that was seen in other peptides
#' @importFrom matrixStats colSums2 colMedians rowMeans2 rowMedians rowSds rowMads
impute_peptide_local_proteinmatrix = function(mat, list_group_indices, mad_prior = NA, par_remove_outliers = TRUE, par_remove_sparse_data = TRUE, return_peptides = FALSE) {
  if(nrow(mat) < 2) {
    return(mat)
  }

  for(g in seq_along(list_group_indices)) {
    g_cols = list_group_indices[[g]]
    index_na = which(is.na(mat[,g_cols,drop=F]), arr.ind = T)
    local_ndatapoints_samples = matrixStats::colSums2(!is.na(mat[,g_cols]))

    # don't use median over samples to scale each peptide, we found that to be quite instable on small datasets, so we use rowmeans instead for finding sample scaling
    # median of means is a pretty robust summary statistic
    local_sample_effect = matrixStats::colMedians(mat[,g_cols] - matrixStats::rowMeans2(mat[,g_cols], na.rm = T), na.rm = T)
    # use the median value for imputing, especially useful for small datasets with outliers
    local_peptide_typical_value = matrixStats::rowMedians(mat[,g_cols], na.rm = T)

    ## optionally, flag extreme values; some value is further from the expected value than N times the estimated peptide's standard deviation
    # NA values are not taken into account, assuming MCAR while finding 'extreme values'
    if(par_remove_outliers) {
      local_peptide_mad = matrixStats::rowSds(mat[,g_cols], na.rm = T)
      # local_peptide_mad = matrixStats::rowMads(mat[,g_cols], na.rm = T)

      local_column_add_sample_effect = local_ndatapoints_samples >= 3 # optionally set to 2
      # iterate local columns
      for(j in which(local_ndatapoints_samples > 0)) {
        # expected value for all peptides; peptide average within-group, and if there are enough data points in column j, take the shift in abundance of these other values into account
        local_col_j__expected_value = local_peptide_typical_value + local_sample_effect[j] * local_column_add_sample_effect[j]

        # compare the 'error' for this peptide to the 'typical MAD' for this protein (use global MAD to smooth out effects if a protein has few peptides and variation is extremely low/high)
        j_error = abs(mat[,g_cols[j]] - local_col_j__expected_value) / mean(c(mad_prior, local_peptide_mad), na.rm = T)
        # alternatively: compare to own MAD
        # j_error = abs(mat[,g_cols[j]] - local_col_j__expected_value) / local_peptide_mad

        # define error threshold
        j_flag = !is.na(j_error) & j_error >= 3
        for(i in which(j_flag)) {
          # for outliers, combine observed value with expected value for outliers. the more extreme the outlier, the less credibility (error used as weight in weighted mean)
          mat[i,g_cols[j]] = weighted.mean(c(mat[i,g_cols[j]], local_col_j__expected_value[i]), w = c(1, j_error[i]))
          # alternatively: simply impute expected value as MCAR, assuming the observed outlier datapoint was a mistake by upstream software
          # mat[i,g_cols[j]] = local_col_j__expected_value[i]
        }
      }
    }

    # finally, impute NAs  (based on values before adjusting extremes)
    if(length(index_na) > 0) {
      for(i in 1:nrow(index_na)) {
        index_peptide = index_na[i,1] # peptide index local/global is the same since the "within group" stuff is just a subset of samples, not a subset of peptides
        index_local_sample = index_na[i,2] # within-group index
        index_global_sample = g_cols[index_local_sample] # index in complete matrix
        # while imputing, always combine peptide's group-average with local sample effect (eg; up/down from typical value in sample s according to shift in observed abundance)
        mat[index_peptide, index_global_sample] = local_peptide_typical_value[index_peptide] + local_sample_effect[index_local_sample]
      }
    }

    ## optionally, drop columns/samples where we have only 1 data point while in other samples we have many (in that sample, imputation heavily relies on a single reference value that is likely an outlier since all other peptides of the same protein are absent)
    # rule: 3+ peptides found in 3+ samples  -->>  any other sample where only 1 peptide was observed is now considered unreliable (eg; don't want to rely on imputation for those samples)
    if(par_remove_sparse_data && sum(local_ndatapoints_samples >= 3) >= 3) {
      mat[,g_cols[local_ndatapoints_samples == 1]] = NA
    }
  }

  if(!return_peptides) {
    mat = log2(matrixStats::colSums2(2^mat, na.rm = T))
  }

  return(mat)
}
