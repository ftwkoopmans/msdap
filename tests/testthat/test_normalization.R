
# validate the normalization function, optimized for speed and making some assumptions while doing so, maintains data integrity


# # very naive log2 data generator
# set.seed(1)
# x = rnorm(100, mean = 5, sd = 1); x[x<1] = 1 # some random data
# m = NULL
# nsamples = 6
# m_groups = letters[1 + (1:nsamples > nsamples/2)] # group labels
# for(i in 1:nsamples) {
#   m = cbind(m, x * rnorm(1, mean = 1, sd = .1) + (i>nsamples/2)) # some noise + bump second group
# }
# m[m<1] = 1 # no negative values, this is supposed to be log2 data
# colnames(m) = 1:ncol(m)
#
# ## add a few percent missing data, MCAR
# m[sample(1:length(m), size = length(m)*0.1)] = NA
#
#
# # plot data
# hist(m)
# boxplot(m)
#
# # to long format for our normalization function
# x = data.table::data.table(m)
# x$key_peptide = 1:nrow(x)
# DT = data.table::melt.data.table(x, id.vars = "key_peptide", variable.name = "key_sample", value.name = "intensity")
# DT$key_sample = as.integer(as.vector(DT$key_sample))
# DT$key_group = m_groups[DT$key_sample] # add groups
# DT[, key_peptide_sample := .GRP, by = list(key_peptide, key_sample)] # generate unique key for peptide*sample
#
# # normalize
# DT$intensity_norm = normalize_intensities(DT, norm_algorithm = "vwmb")
#
# # qc plots
# DT_wide_prenorm = data.table::dcast(DT, key_peptide ~ key_sample, value.var = "intensity", fill = NA)
# DT_wide_postnorm = data.table::dcast(DT, key_peptide ~ key_sample, value.var = "intensity_norm", fill = NA)
#
# # input matrix is the same as after casting to long, then to wide again
# all(m == as.matrix(DT_wide_prenorm)[,-1], na.rm = T)
# # boxplot pre-normalization (same as boxplot(m) above)
# boxplot(as.matrix(DT_wide_prenorm)[,-1])
# # after normalization, all levels should be similar
# boxplot(as.matrix(DT_wide_postnorm)[,-1])
