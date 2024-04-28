
#' check if peptide tibble has cache
#' @param dataset your dataset
#' @export
check_dataset_hascache = function(dataset) {
  is.list(dataset) && all(c("groups","dt_pep_group") %in% names(dataset))
}



#' invalidate peptide tibble cache
#' @param dataset your dataset
#' @export
invalidate_cache = function(dataset) {
  dataset$groups = NULL
  dataset$dt_pep_group = NULL
  return(dataset)
}



#' Cache data required for fast downstream peptide filtering
#'
#' eg; number of replicate samples per group in which a peptide is detected, within-group variation, etc.
#'
#' @param dataset your dataset
#' @importFrom data.table data.table setDT setkey chmatch merge.data.table melt.data.table
#' @export
cache_filtering_data = function(dataset) {
  start_time = Sys.time()
  check_dataset_integrity(dataset)

  # add keys to samples table
  s = data.table::setDT(dataset$samples) # convert to data.table, editing this object BY REFERENCE downstream
  s[, key_sample := .GRP, by = sample_id]
  s[, key_group := .GRP, by = group]
  rm(s)

  # add keys to peptides table
  i = data.table::chmatch(dataset$peptides$sample_id, dataset$samples$sample_id)
  dataset$peptides$key_sample = dataset$samples$key_sample[i]
  dataset$peptides$key_group = dataset$samples$key_group[i]
  rm(i)

  # new reference table for groups
  grps = dataset$samples %>% group_by(group) %>% summarise(size=n(), size_noexclude=sum(!exclude), key_group=key_group[1]) %>% arrange(key_group)

  # use a subset of peptides where need to make unique keys
  p = data.table::setDT(dataset$peptides) # convert to data.table, editing this object BY REFERENCE downstream
  p[, key_peptide := .GRP, by = peptide_id]
  p[, key_protein := .GRP, by = protein_id]
  p[, key_peptide_sample := .GRP, by = list(key_peptide, key_sample)]
  p[, key_protein_group := .GRP, by = list(key_protein, key_group)]
  p[, key_peptide_group := .GRP, by = list(key_peptide, key_group)]
  # flag 'exclude' samples
  kpg = p$key_peptide_group
  kpg[p$key_sample %in% (dataset$samples %>% filter(exclude==T) %>% pull(key_sample))] = NA
  p[, key_peptide_group_noexclude := kpg]
  # detect not as a boolean, but scaled by number of detect per sample
  # ` * !is.na(key_peptide_group_noexclude)` is simply to remove all data points for samples that are excluded, which allows us to downstream group by key_group and not worry about excluded samples
  p[, `:=` (detect_scaled_noexclude = (detect/sum(detect)) * !is.na(key_peptide_group_noexclude),
            quant_scaled_noexclude = 1/.N * !is.na(key_peptide_group_noexclude)),
    by = key_sample]


  ## efficiently collapse filtering data for each peptide * 'sample group'
  dt_pep_group = p[ , .(key_peptide=key_peptide[1], key_protein=key_protein[1], key_group=key_group[1], key_protein_group=key_protein_group[1],
                        nquant = .N,
                        ndetect = sum(detect),
                        nquant_noexclude = sum(!is.na(key_peptide_group_noexclude)),
                        ndetect_noexclude = sum(!is.na(key_peptide_group_noexclude) & detect),
                        ndetect_scaled_noexclude = sum(detect_scaled_noexclude),
                        nquant_scaled_noexclude = sum(quant_scaled_noexclude)),
                    by=key_peptide_group]
  rm(p)

  # add group sizes and cache the 'fraction of samples within each group' where a peptide is detected
  data.table::setkey(dt_pep_group, key_group) # set key for data.table::merge()
  dt_criteria = data.table::data.table(grps %>% select(key_group, size, size_noexclude) %>% mutate_all(as.integer), key="key_group")
  dt_pep_group = data.table::merge.data.table(dt_pep_group, dt_criteria, all.x = T, sort = FALSE)
  dt_pep_group[ , `:=`(ndetect_fraction = ndetect / size, ndetect_noexclude_fraction = ndetect_noexclude / size_noexclude)]


  ################# pre-cache CoV
  ###### simplest approach: just compute CoV for each peptide, and later worry about the subset we will use downstream yes/no. These are only reference CoV's for sorting peptides @ topN filter anyway
  m = as_matrix_except_first_column(dataset$peptides %>% pivot_wider(id_cols = "key_peptide", names_from = "key_sample", values_from = "intensity"))
  # remove rows that lack data
  m = m[matrixStats::rowSums2(!is.na(m)) >= 2, , drop=F]
  if(nrow(m) > 1) {
    # fast mode normalization. samples$group is a character array, enforced upstream by sample_metadata_sort_and_filter()
    m = normalize_matrix(m, group_by_cols = dataset$samples$group[match(colnames(m), as.character(dataset$samples$key_sample))], algorithm = "vwmb")
    # debug; boxplot(m, outline=F)
    # convert log2 intensities to natural log for CoV
    m = log2_to_ln(m)

    # cache CoV (so we get accurate estimation of 'low variation peptides' in downstream topN filters)
    dtw_cov = data.table::data.table(key_peptide = as.integer(rownames(m)))
    for(k in grps$key_group) {
      sid = intersect(as.character(dataset$samples %>% filter(key_group == k & !exclude) %>% pull(key_sample)), colnames(m))
      if(length(sid) >= 2) {
        dtw_cov[,as.character(k)] = coefficient_of_variation_vectorized(m[,colnames(m) %in% sid,drop=F])
      }
    }
    dt_cov = data.table::melt.data.table(dtw_cov, id.vars = "key_peptide", variable.name = "key_group", value.name = "cov")
  }

  ## convert CoV estimates to scores between 0~1
  # this parameter has huge impact on peptide topN selection and is critical.
  #
  # Using the ecdf() seems like a natural choice, but on tiny datasets with lots of ties this leads to problems, such as; ecdf(c(1,1,1,2,3))(c(1,1,1,2,3)); ecdf(c(1,2,3))(c(1,2,3))
  # instead, threshold the top and bottom 1% quantiles, then scale between 0~1.
  # CoV's are log-normal. we could relax penalizing 'outlier' CoVs by log transforming the CoVs before scaling 0~1
  #
  # alternatively, we could scale between 0.1~min(1, quantile(dt_cov$cov,0.99,na.rm=T)) for some added robustness. Eg; if a dataset has barely any variation, a peptides with ~5% variation is 'much worse' than one at 2% while at that point, difference in #detect might weigh more more. Haven't found any such problems in real data though
  if(exists("dt_cov") && "cov" %in% colnames(dt_cov)) {
    dt_cov$cov_scale_quantiles = scale_between_quantiles(dt_cov$cov, min_quantile = 0.01, max_quantile = 0.99)
    dt_cov$key_group = as.integer(as.vector(dt_cov$key_group))
    # join CoV's into peptide*group tibble
    data.table::setkey(dt_cov, key_peptide, key_group)
    data.table::setkey(dt_pep_group, key_peptide, key_group)
    dt_pep_group = data.table::merge.data.table(dt_pep_group, dt_cov, all.x = T, sort = FALSE)
  } else {
    dt_pep_group[, cov_scale_quantiles := NA ]
  }

  ## peptide*group 'quality scores'
  # the 'exclude' samples are completely disregarded while computing this score
  # combined score, where values missing from either score are set to zero, used downstream to select topN peptides; fraction detect + 1-cov_score. cov scores are inversed, because lower CoV = better score
  dt_pep_group[, score_detect_cov := replace_nonfinite(ndetect_noexclude_fraction, 0) + replace_nonfinite(1 - cov_scale_quantiles, 0) ]

  dataset$dt_pep_group = dt_pep_group
  dataset$groups = grps
  # guarantee downstream compatability, enforce tibble type
  dataset$peptides = as_tibble(dataset$peptides)
  dataset$samples = as_tibble(dataset$samples)
  dataset$proteins = as_tibble(dataset$proteins)
  append_log_timestamp("caching filter data", start_time)
  return(dataset)

  ########################################################################################################################################
  #
  # mat_group_cov = matrix(NA, nrow = nrow(m), ncol=nrow(grps), dimnames=list(rownames(m), grps$key_group))
  # for(k in grps$key_group) {
  #   sid = intersect(as.character(dataset$samples %>% filter(key_group == k & !exclude) %>% pull(key_sample)), colnames(m))
  #   if(length(sid) >= 2) {
  #     mat_group_cov[,as.character(k)] = coefficient_of_variation_vectorized(m[,colnames(m) %in% sid,drop=F])
  #   }
  # }
  #debug; boxplot(mat_group_cov, names = grps$group[match(colnames(mat_group_cov), as.character(grps$key_group))], ylim=c(0,2), outline=F, las=2)
  #
  # #cleanup
  # rm(dt_cov)
  # rm(dtw_cov)
  # rm(m)
  #
  ###### backup code; select only peptides with N detect per group for CoV computation
  # # filter rules; which peptides are we computing CoV for?
  # dt_filter_group[ , `:=`(pass = ndetect>=1 & nquant>=2, pass_noexclude = ndetect_noexclude>=1 & nquant_noexclude>=2)]
  #
  # ## non-contrast filters
  # # note; never ifelse large arrays, copy entire array then mask to replace values with NA (see benchmark code below)
  # dataset$peptides$intensity_by_group_temp = dataset$peptides$intensity
  # # add 'by group' filter ignoring the exclude samples where we remove from the intensity values those that originate from a sample where the respective peptide doesn't pass group-wide filter (no normalization applied yet)
  # i = match(dataset$peptides$key_peptide_group, dt_filter_group$key_peptide_group)
  # dataset$peptides$intensity_by_group_temp[is.na(i) | !dt_filter_group$pass_noexclude[i]] <- NA
  # # cleanup
  # rm(i)
  # rm(dt_filter_group)
  #
  # #### we can pre-cache CoV per peptide. Only difference with computing during filter is slight influence on normalization, but CoV estimate is impacted for all peptides anyway. relative CoV ranking is very similar pre/post filter (for the same peptide)
  # ## cache data for topN filter; within-group variation  =  after by-group filter + norm, do fast CoV
  # # prep data in wide format
  # dt_int = data.table::data.table(dataset$peptides %>%
  #                                   select(key_peptide, key_group, key_sample, intensity=intensity_by_group_temp) %>%
  #                                   filter(!is.na(intensity)) %>%
  #                                   # !! non-log intensities for norm
  #                                   mutate(intensity = 2^intensity))
  # dtw_int = data.table::dcast(dt_int, key_peptide ~ key_sample, value.var = "intensity", fill = NA)
  # # fast mode normalization
  # m = as.matrix(dtw_int[,-"key_peptide"])
  # rownames(m) = dtw_int$key_peptide
  # m = log2_to_ln(normalize_mode(m, group_by_cols = dataset$samples$group[match(colnames(m), as.character(dataset$samples$key_sample))]))
  #
  # #cleanup
  # rm(dtw_int)
  # rm(m)
  # peptides$intensity_by_group_temp = NULL
}



#' Filter dataset
#'
#' For each of the filters (by group, all groups, by contrast) an extra column will be appended to the dataset$peptides table that contains the intensity data for all peptides that pass these filters (`name: intensity_<filter>`).
#' Optionally, you can apply normalization (recommended).
#'
#' Note; this is built-in for \code{analysis_quickstart} so if you use that all-in-one function you don't have to perform all pipeline steps manually
#'
#' @param dataset a valid dataset object generated upstream by, for instance, import_dataset_skyline
#' @param filter_min_detect in order for a peptide to 'pass' in a sample group, in how many replicates must it be detected?
#' @param filter_fraction_detect in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be detected?
#' @param filter_min_quant in order for a peptide to 'pass' in a sample group, in how many replicates must it be quantified?
#' @param filter_fraction_quant in order for a peptide to 'pass' in a sample group, what fraction of replicates must it be quantified?
#' @param filter_min_peptide_per_prot in order for a peptide to 'pass' in a sample group, how many peptides should be available after detect filters?
#' @param filter_topn_peptides maximum number of peptides to maintain for each protein (from the subset that passes above filters, peptides are ranked by the number of samples where detected and their variation between replicates). 1 is default, 2 can be a good choice situationally. If set to 1, make sure to inspect individual peptide data/plots for proteins with 1 peptide.
#' @param norm_algorithm normalization algorithms. options; "", "vsn", "loess", "rlr", "msempire", "vwmb", "modebetween". Provide an array of options to run each algorithm consecutively
#' @param rollup_algorithm the algorithm for combining peptides to proteins as used in normalization algorithms that require a priori rollup from peptides to a protein-level abundance matrix (e.g. modebetween_protein). Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#' @param by_group within each sample group, apply the filter. All peptides that fail the filter in group g will have intensity value NA in the intensity_by_group column for the samples in the respective group
#' @param all_group in every sample group, apply the filter. All peptides that fail the filter in any group will have intensity value NA in the intensity_all_groups column for all samples
#' @param by_contrast should the above filters be applied to all sample groups, or only those tested within each contrast? Enabling this optimizes available data in each contrast, but increases the complexity somewhat as different subsets of peptides are used in each contrast and normalization is applied separately
#'
#' @importFrom data.table data.table setDT setkey merge.data.table dcast
#' @importFrom matrixStats rowSums2
#' @export
filter_dataset = function(dataset,
                          # peptide filter criteria applied within each sample group
                          filter_min_detect = 1, filter_fraction_detect = 0, filter_min_quant = 0, filter_fraction_quant = 0,
                          # respective criteria on protein level
                          filter_min_peptide_per_prot = 1, filter_topn_peptides = 0,
                          # normalization
                          norm_algorithm = "",
                          rollup_algorithm = "maxlfq",
                          # which filters to apply
                          by_group = T, all_group = T, by_contrast = F) {

  start_time = Sys.time()

  ##### input validation
  # validate function parameters
  check_parameter_is_numeric(filter_min_detect, filter_fraction_detect, filter_min_quant, filter_fraction_quant, filter_min_peptide_per_prot, filter_topn_peptides)
  check_parameter_is_boolean(by_group, all_group, by_contrast)
  if(!all(is.character(norm_algorithm))) {
    append_log("function argument 'norm_algorithm' must be a string", type = "error")
  }
  if(length(rollup_algorithm) != 1 || any(!is.character(rollup_algorithm))) {
    append_log("function argument 'rollup_algorithm' must be a single string (not an array)", type = "error")
  }

  # there should be no decoys at this point
  if("isdecoy" %in% colnames(dataset$peptides) && any(dataset$peptides$isdecoy)) {
    append_log("decoys should be removed prior to downstream data analysis. Either exclude them while importing using import_dataset_x() as is default, or remove them from the dataset prior to calling this function. eg; dataset$peptides = dataset$peptides %>% filter(isdecoy == FALSE)", type = "error")
  }

  # there should be no sample fractions at this point
  if ("fraction" %in% colnames(dataset$samples)) {
    append_log("fractionated data is provided to the filter_dataset() function, but fractions should have been merged prior! You probably have to run 'dataset = merge_fractionated_samples(dataset)', please refer to the online documentation or check the implementation of the analysis_quickstart() function for example code.", type = "error")
  }

  check_dataset_integrity(dataset)
  ##### input validation


  if(!check_dataset_hascache(dataset)) {
    dataset = cache_filtering_data(dataset)
  }

  # used for reporting results later on
  npep_input = n_distinct(dataset$peptides$peptide_id)
  log_stats = NULL

  # remove all pre-existing filtering columns
  cols_intensity = grep("^intensity_", colnames(dataset$peptides), ignore.case = T, value=T)
  if(length(cols_intensity) > 0) {
    append_log(sprintf("applying new filters, removing pre-existing data columns: %s", paste(cols_intensity, collapse=", ")))
    dataset$peptides[,cols_intensity] = NULL # remove additional intensity columns
    dataset$peptides = dataset$peptides %>% filter(is.finite(intensity)) # remove imputed values, if any
  }


  any_samples_excluded = any(dataset$samples$exclude)
  ## translate filters to criteria per group; 'detect in x% of samples' filter depend on group sizes
  dataset$groups$ndetect_min = pmin(dataset$groups$size, pmax(filter_min_detect, ceiling(dataset$groups$size * filter_fraction_detect)))
  dataset$groups$nquant_min  = pmin(dataset$groups$size, pmax(filter_min_quant,  ceiling(dataset$groups$size * filter_fraction_quant)))
  dataset$groups$ndetect_noexclude_min = pmin(dataset$groups$size_noexclude, pmax(filter_min_detect, ceiling(dataset$groups$size_noexclude * filter_fraction_detect)))
  dataset$groups$nquant_noexclude_min  = pmin(dataset$groups$size_noexclude, pmax(filter_min_quant,  ceiling(dataset$groups$size_noexclude * filter_fraction_quant)))

  ## apply filter to each group
  # first, join the requirements
  data.table::setkey(dataset$dt_pep_group, key_group) # set key for data.table::merge()
  dt_criteria = data.table::data.table(dataset$groups %>% select(key_group, ndetect_min, nquant_min, ndetect_noexclude_min, nquant_noexclude_min) %>% mutate_all(as.integer), key="key_group")
  dt_filter_group = data.table::merge.data.table(dataset$dt_pep_group, dt_criteria, all.x = T, sort = FALSE)
  dt_filter_group[ , `:=`(pass = ndetect >= ndetect_min & nquant >= nquant_min,
                          pass_noexclude = ndetect_noexclude >= ndetect_noexclude_min & nquant_noexclude >= nquant_noexclude_min)]

  ## now we know for each peptide whether it passes the filter criteria in each group (with and without taking 'exclude' samples into account)

  ## wide format for re-use (eg; by-group or all-group)
  data.table::setkey(dt_filter_group, pass_noexclude) # set key for filter step below; remove all peptide*group that don't pass anywhere to minimize wide-table size
  dt_filter_group_wide_noexclude = data.table::dcast(dt_filter_group[pass_noexclude>0], key_peptide ~ key_group, value.var = "pass_noexclude", fill = FALSE)
  dt_filter_group_wide_noexclude[ , `:=`(ngroup = matrixStats::rowSums2(as.matrix(.SD))), .SDcols=eval(quote(setdiff(colnames(dt_filter_group_wide_noexclude), "key_peptide")))]
  # analogous
  data.table::setkey(dt_filter_group, pass)
  dt_filter_group_wide = data.table::dcast(dt_filter_group[pass>0], key_peptide ~ key_group, value.var = "pass", fill = FALSE)
  dt_filter_group_wide[ , `:=`(ngroup = matrixStats::rowSums2(as.matrix(.SD))), .SDcols=eval(quote(setdiff(colnames(dt_filter_group_wide), "key_peptide")))]


  if(by_group) {
    # the by-group filter is a special-case; separately filter each group (detect, minpep, topN), then combine everything in 1 column (by definition, there is no overlap because each sample belongs to exactly 1 sample group)
    i = match(dataset$peptides$key_peptide_group, dt_filter_group$key_peptide_group)
    for(k in dataset$groups$key_group) {
      n = paste0("intensity_by_group_", k)
      dataset$peptides[,n] = dataset$peptides$intensity
      # subset of rows matching current group
      rows_k = dataset$peptides$key_group == k
      # a) we don't want intensity values for 'exclude' samples in the result. b) use `dt_filter_group$pass_noexclude` to check if a peptide passes filters in each peptide*samplegroup_k
      dataset$peptides[!rows_k | is.na(dataset$peptides$key_peptide_group_noexclude) | is.na(i) | !dt_filter_group$pass_noexclude[i], n] <- NA
      # debug; print(n);print(dataset$peptides %>% select(value = !!n, key_group, sample_id, peptide_id) %>% filter(!is.na(value)) %>% count(sample_id) %>% arrange(n) %>% left_join(dataset$samples %>% select(sample_id, group, exclude)), n=99)    }
    }
    rm(i)
  }


  if(all_group) {
    ## global filter; only peptides that pass each group
    # note; never ifelse large arrays, copy entire array then mask to replace values with NA (see benchmark code below)
    dataset$peptides$intensity_all_group = dataset$peptides$intensity
    # a) we don't want intensity values for 'exclude' samples in the result. b) use `dt_filter_group_wide$ngroup` to check if a peptide passes filters in all tested peptide*samplegroup
    i = match(dataset$peptides$key_peptide, dt_filter_group_wide_noexclude$key_peptide)
    dataset$peptides$intensity_all_group[is.na(dataset$peptides$key_peptide_group_noexclude) | is.na(i) | dt_filter_group_wide_noexclude$ngroup[i] != nrow(dataset$groups)] <- NA
    # analogous
    if(any_samples_excluded) {
      dataset$peptides$intensity_all_group_withexclude = dataset$peptides$intensity
      i = match(dataset$peptides$key_peptide, dt_filter_group_wide$key_peptide)
      dataset$peptides$intensity_all_group_withexclude[is.na(i) | dt_filter_group_wide$ngroup[i] != nrow(dataset$groups)] <- NA
      rm(i)
    }
    # log results
    log_stats = c(log_stats, sprintf("%d/%d peptides were retained after filtering over all groups", dataset$peptides %>% filter(!is.na(intensity_all_group)) %>% distinct(peptide_id) %>% nrow(), npep_input))
  }


  #### filter and normalize by contrast
  # get all contrasts from sample table
  column_contrasts = dataset_contrasts(dataset)
  # if there are no contrasts setup, disable 'by_contrast'
  by_contrast = by_contrast && length(column_contrasts) > 0

  if(by_contrast) {
    for (col_contr in column_contrasts) { # col_contr = column_contrasts[1]
      # only keep the current contrast in samples table, so downstream code tests this contrast only
      # note; outlier samples are not included in contrasts as at all in upstream 'contrasts implementation'
      col_contr_samples = dataset$samples %>% select(key_sample, key_group, group, contrast = !!col_contr) %>% filter(contrast != 0)
      grp1_key = col_contr_samples %>% filter(contrast == 1) %>% distinct(key_group) %>% pull()
      grp2_key = col_contr_samples %>% filter(contrast == 2) %>% distinct(key_group) %>% pull()
      # grp1 = col_contr_samples %>% filter(contrast == 1) %>% pull(group)
      # grp2 = col_contr_samples %>% filter(contrast == 2) %>% pull(group)
      # grp1_key = dataset$groups %>% filter(group %in% grp1) %>% distinct(key_group) %>% pull()
      # grp2_key = dataset$groups %>% filter(group %in% grp2) %>% distinct(key_group) %>% pull()

      if(length(grp1_key) == 1) {
        grp1_mask = dataset$peptides$key_peptide %in% dt_filter_group_wide_noexclude$key_peptide[ dt_filter_group_wide_noexclude[[as.character(grp1_key)]] ]
      } else {
        # if multi-group, select all columns for this group and test that all are true
        dt_filter_group_wide_noexclude[, `:=`(temp = matrixStats::rowSums2(as.matrix(.SD))), .SDcols=eval(quote(as.character(grp1_key)))]
        grp1_mask = dataset$peptides$key_peptide %in% dt_filter_group_wide_noexclude$key_peptide[ dt_filter_group_wide_noexclude$temp == length(grp1_key) ]
      }

      if(length(grp2_key) == 1) {
        grp2_mask = dataset$peptides$key_peptide %in% dt_filter_group_wide_noexclude$key_peptide[ dt_filter_group_wide_noexclude[[as.character(grp2_key)]] ]
      } else {
        # if multi-group, select all columns for this group and test that all are true
        dt_filter_group_wide_noexclude[, `:=`(temp = matrixStats::rowSums2(as.matrix(.SD))), .SDcols=eval(quote(as.character(grp2_key)))]
        grp2_mask = dataset$peptides$key_peptide %in% dt_filter_group_wide_noexclude$key_peptide[ dt_filter_group_wide_noexclude$temp == length(grp2_key) ]
      }


      intensity_col_contr = paste0("intensity_", col_contr)
      dataset$peptides[,intensity_col_contr] = dataset$peptides$intensity
      dataset$peptides[!(grp1_mask & grp2_mask & dataset$peptides$key_sample %in% col_contr_samples$key_sample), intensity_col_contr] <- NA

      # log results
      log_stats = c(log_stats, sprintf("%d/%d peptides were retained after filtering within %s", dataset$peptides %>% select(peptide_id, intensity = !!intensity_col_contr) %>% filter(!is.na(intensity)) %>% distinct(peptide_id) %>% nrow(), npep_input, col_contr))

      # table(grp1_mask); table(grp2_mask); table(!(grp1_mask & grp2_mask & dataset$peptides$key_group %in% c(grp1_key,grp2_key))); dataset$peptides %>% select(value = !!intensity_col_contr, key_group, sample_id, peptide_id) %>% filter(!is.na(value)) %>% count(key_group, sample_id)
    }
  # } else {
  #   dataset$peptides[,paste0("intensity_", column_contrasts)] = NULL # drop column if pre-existing. eg; suppose user changes filters + normalization, then doesn't want this filter; cannot keep old data
  }



  ###### auto-apply minpep- and topn-filters to all 'intensity_' columns and finally normalise

  cols_intensity = grep("^intensity_", colnames(dataset$peptides), ignore.case = T, value=T)

  # apply min-peptide-count filter to all filter columns
  if(filter_min_peptide_per_prot > 1) {
    for(col_intensity in cols_intensity) { # col_intensity=cols_intensity[2]
      key_protein_valid = filter_proteins_by_pepcount(dataset$peptides, col_intensity=col_intensity, min_peptides = filter_min_peptide_per_prot)
      dataset$peptides[!(dataset$peptides$key_protein %in% key_protein_valid), col_intensity] <- NA # invalidate those peptides that are not in the 'valid protein set'
    }
  }


  # apply topn_peptides filter to all filter columns
  if(filter_topn_peptides > 0) {
    for(col_intensity in cols_intensity) {
      key_peptide_valid = filter_proteins_by_topn(dataset$peptides, dt_pep_group=dataset$dt_pep_group, col_intensity=col_intensity, filter_topn_peptides = filter_topn_peptides)
      dataset$peptides[!(dataset$peptides$key_peptide %in% key_peptide_valid), col_intensity] <- NA # invalidate those peptides that are not in the 'valid peptide set'
    }
  }


  # merge all separate by-group filters into one column
  if(by_group) {
    cols_bygroup = grep("^intensity_by_group_", colnames(dataset$peptides), value=T)
    dataset$peptides$intensity_by_group = NA
    for(k in cols_bygroup) {
      dataset$peptides$intensity_by_group = pmin(dataset$peptides %>% pull(!!k), dataset$peptides$intensity_by_group, na.rm = TRUE) # leaves tibble intact, unlike mutate() which resets tibble attributes
      # print(sum(!is.na(dataset$peptides$intensity_by_group))) #debug
    }
    dataset$peptides = dataset$peptides %>% select(-one_of(cols_bygroup))
    # update intensity column names after merging by-group filter
    cols_intensity = grep("^intensity_", colnames(dataset$peptides), ignore.case = T, value=T)
    # log results
    log_stats = c(log_stats, sprintf("%d/%d peptides were retained after filtering within each group independently (\"by group\")", dataset$peptides %>% filter(!is.na(intensity_by_group)) %>% distinct(peptide_id) %>% nrow(), npep_input))
  }


  if(any(norm_algorithm != "")) {
    for(col_intensity in cols_intensity) { #col_intensity=cols_intensity[1]
      dataset = normalize_peptide_intensity_column(dataset, col_intensity = col_intensity, norm_algorithm = norm_algorithm, rollup_algorithm = rollup_algorithm)
      # dataset$peptides[ , col_intensity] = normalize_intensities(data.table::setDT(dataset$peptides %>% select(key_peptide_sample, key_peptide, key_protein, key_sample, key_group, intensity = !!col_intensity)), norm_algorithm = norm_algorithm)
    }
  }


  # report to console
  if(length(log_stats) > 0) {
    log_settings = NULL
    if(filter_min_detect > 0) log_settings = c(log_settings, paste("min_detect =", filter_min_detect))
    if(filter_fraction_detect > 0) log_settings = c(log_settings, paste("fraction_detect =", filter_fraction_detect))
    if(filter_min_quant > 0) log_settings = c(log_settings, paste("min_quant =", filter_min_quant))
    if(filter_fraction_quant > 0) log_settings = c(log_settings, paste("fraction_quant =", filter_fraction_quant))
    if(filter_min_peptide_per_prot > 1) log_settings = c(log_settings, paste("min_peptide_per_prot =", filter_min_peptide_per_prot))
    if(filter_topn_peptides > 0) log_settings = c(log_settings, paste("topn_peptides =", filter_topn_peptides))
    if(any(norm_algorithm != "")) log_settings = c(log_settings, paste0("norm_algorithm = '", paste(norm_algorithm, collapse = "&"), "'"))
    log_settings = c(log_settings, paste0("rollup_algorithm = '", rollup_algorithm, "'"))

    append_log(sprintf("filter dataset with settings: %s\n%s", paste(log_settings, collapse = "; "), paste(log_stats, collapse = "\n")), type = "info")
  }

  # guarantee downstream compatability, enforce tibble type
  dataset$peptides = as_tibble(dataset$peptides)
  dataset$samples = as_tibble(dataset$samples)
  dataset$proteins = as_tibble(dataset$proteins)
  append_log_timestamp("peptide filtering and normalization", start_time)
  return(dataset)

  #################  benchmark code; don't ifelse()
  # rbenchmark::benchmark(a={dataset$peptides$intensity_contrast_test = ifelse(dataset$peptides$key_peptide %in% dt_filter_group_wide$key_peptide[grp1_mask & grp2_mask], dataset$peptides$intensity, NA)},
  #                       b={dataset$peptides$intensity_contrast_test = dataset$peptides$intensity; dataset$peptides$intensity_contrast_test[grp1_mask & grp2_mask] <- NA},
  #                       replications = 100)
  # result; a=12sec, b=1sec
  #
  ################# some reference code
  ## ref; all-in-one bygroup filter & minpep
  # dataset$peptides$intensity_by_group = dataset$peptides$intensity
  # # add 'by group' filter ignoring the exclude samples where we remove from the intensity values those that originate from a sample where the respective peptide doesn't pass group-wide filter (no normalization applied yet)
  # i = match(dataset$peptides$key_peptide_group, dt_filter_group$key_peptide_group)
  # # a) we don't want intensity values for 'exclude' samples in the result. b) use `dt_filter_group$pass_noexclude` to check if a peptide passes filters in each peptide*samplegroup
  # dataset$peptides$intensity_by_group[is.na(dataset$peptides$key_peptide_group_noexclude) | is.na(i) | !dt_filter_group$pass_noexclude[i]] <- NA
  # if(filter_min_peptide_per_prot > 1) {
  #   # need an extra level of summary to count, within each group, how many peptides pass the filter; sum(pass_noexclude)>=N   by {group,protein}
  #   tmp = dt_filter_group[ , .(pass=sum(pass_noexclude)>=filter_min_peptide_per_prot), by=key_protein_group][pass==T]
  #   dataset$peptides$intensity_by_group[!(dataset$peptides$key_protein_group %in% tmp$key_protein_group)] <- NA
  #   rm(tmp)
  # }
  # rm(i)
  #
  #
  # # 'binary encoding' to combine both filter outcomes in one variable (so we don't need multiple wide-format matrices downstream)
  # # example code; for(a in c(F,T)) { for(b in c(F,T)) { print(c(as.character(a), as.character(b), a+b*2)) }}
  # # reference; 0 = doesn't pass anywhere. 1=pass&!pass_noexclude 2=!pass&pass_noexclude 3=pass&pass_noexclude
  # dt_filter_group[ , `:=`(pass_code = pass + pass_noexclude*2)]
  # ## wide format for re-use (eg; by-group or all-group)
  # data.table::setkey(dt_filter_group, pass_code) # set key for filter step below; remove all peptide*group that don't pass anywhere to minimize wide-table size
  # dt_filter_group_wide = data.table::dcast(dt_filter_group[pass_code>0], key_peptide ~ key_group, value.var = "pass_code", fill = 0)
}







###########




filter_proteins_by_pepcount = function(peptides, col_intensity, min_peptides) {
  protein_filter = peptides %>%
    select(key_peptide, key_protein, value=!!col_intensity) %>%
    filter(!is.na(value)) %>%
    select(key_peptide, key_protein) %>%
    distinct(key_peptide, .keep_all = T) %>%
    count(key_protein, name = "peptide_count")

  # log_count_pre = nrow(protein_filter)
  protein_filter = protein_filter %>%
    filter(peptide_count >= min_peptides)
  # log_count_post = nrow(protein_filter)

  return(protein_filter$key_protein)
  # peptides[!(peptides$key_protein %in% protein_filter$key_protein), col_intensity] <- NA

  ## reference code for both dplyr and data.table in a benchmark. Performance is comparable within a few percent
  # col_intensity = "intensity_all_group"
  # col_intensity = "intensity"
  # min_peptides=2
  # rbenchmark::benchmark(datatable_1 = {
  #   x = unique(data.table::data.table(peptides %>% select(key_peptide, key_protein, value=!!col_intensity))[!is.na(value)], by = "key_peptide")[ , peptide_count := .N, by=key_protein][peptide_count >= min_peptides]
  # },
  # datatable_2 = {
  #   x = unique(data.table::data.table(peptides %>% select(key_peptide, key_protein, value=!!col_intensity))[!is.na(value)], by = "key_peptide")[ , .( pass= .N >= min_peptides), by=key_protein][pass==T]
  # },
  # datatable_3 = {
  #   x = unique(data.table::data.table(peptides %>% select(key_peptide, key_protein, value=!!col_intensity) %>% filter(!is.na(value)) %>% select(-value)), by = "key_peptide")[ , peptide_count := .N, by=key_protein][peptide_count >= min_peptides]
  # },
  # dplyr = {
  #   protein_filter = peptides %>%
  #     select(key_peptide, key_protein, value=!!col_intensity) %>%
  #     filter(!is.na(value)) %>%
  #     select(key_peptide, key_protein) %>%
  #     distinct(key_peptide, .keep_all = T) %>%
  #     count(key_protein, name = "peptide_count") %>%
  #     filter(peptide_count >= min_peptides)
  # },
  # replications = 100)
  # double-check both solutions provide same results
  # all(x$key_protein %in% protein_filter$key_protein)
  # all(protein_filter$key_protein %in% x$key_protein)
}



##
# 1) given the subset of peptides that survived filters so far
# 2) select topN best per protein based on #detect (larger is better) and CoV (lower is better)
## filter:
# 1) subset peptides that pass previous filters
# 2) sort by combined score (sort decreasing=T)
# 3) subset topN
# 4) return valid peptides
#' @importFrom data.table setDT setkey merge.data.table setorder
filter_proteins_by_topn = function(peptides, dt_pep_group, col_intensity, filter_topn_peptides = 10) {
  # long-format table of peptide*group keys, only for those peptides that passes earlier filters (!!col_intensity still has a value)
  # !! make sure peptides are not doubly counted -->> take distinct `key_peptide_group`

  # now we have a unique set of peptide*group that pass filters @ col_intensity
  x = unique(data.table::setDT(peptides %>% select(key_peptide_group, key_peptide, key_protein, intensity=!!col_intensity) %>%
                                      # only use the intensity column to select peptide*sample that have not been previously removed
                                      filter(!is.na(intensity)) %>% select(-intensity), key = "key_peptide_group"),
             by = "key_peptide_group")

  # merge only the `score_detect_cov` from `dt_pep_group` in a left-join
  data.table::setkey(dt_pep_group, key_peptide_group)
  x = data.table::merge.data.table(x, dt_pep_group[ , .(key_peptide_group, score_detect_cov)], all.x = T, sort = FALSE)

  # summarise scores per peptide: sum scores over all unique key_peptide_group
  x = x[ , .(key_protein=key_protein[1], score_detect_cov=sum(score_detect_cov, na.rm=T)), by=key_peptide]

  # sort by score, then take topN
  data.table::setorder(x, -score_detect_cov, na.last = T) # example code to double-check sorting: x = data.table(a=1:3, b=c(F,T,F));x;setorder(x, -b,-a, na.last=FALSE);x
  # 1) number peptides within a protein from 1:N, 2) remove those outside of topN as an efficient way of subsetting
  x = x[ , index := seq_len(.N), by = key_protein][index <= filter_topn_peptides]

  return(x$key_peptide)

  ## reference code for both dplyr and data.table in a benchmark. Performance is comparable, with a slight edge for data.table
  # rbenchmark::benchmark(dplyr={
  #   x = data.table::data.table(peptides %>%
  #                                select(key_peptide_group, key_peptide, key_protein, intensity=!!col_intensity) %>%
  #                                # only use the intensity column to select peptide*sample that have not been previously removed
  #                                filter(!is.na(intensity)) %>%
  #                                select(-intensity) %>%
  #                                # now we have a unique set of peptide*group that pass filters @ col_intensity
  #                                distinct(key_peptide_group, .keep_all = T),
  #                              key="key_peptide_group")
  # },
  # datatable_1 = { x = unique(data.table::data.table(peptides %>% select(key_peptide_group, key_peptide, key_protein, intensity=!!col_intensity) %>% filter(!is.na(intensity)) %>% select(-intensity)), by = "key_peptide_group") },
  # datatable_2 = { x = unique(data.table::data.table(peptides %>% select(key_peptide_group, key_peptide, key_protein, intensity=!!col_intensity))[!is.na(intensity)], by = "key_peptide_group") },
  # replications = 100
  # )
}



# helper function for simple filtering, applies a filter to each group individually
# eg; within each group, remove peptides that are not detected in at least 2 replicates
#' @importFrom data.table data.table
filter_peptides_by_group = function(dataset, colname="intensity_temp", disregard_exclude_samples=F, nquant=2, ndetect=1, norm_algorithm = "vwmb", rollup_algorithm = "maxlfq") {
  # input validation
  check_dataset_integrity(dataset)
  if(length(colname) != 1 || !grepl("^intensity_", colname)) {
    append_log(paste("the column name that should hold resulting peptide intensities must start with 'intensity_' to ensure compatability with the rest of the codebase. Provided value:", colname), type = "error")
  }
  if(!check_dataset_hascache(dataset)) {
    dataset = cache_filtering_data(dataset)
  }

  # align with pre-cached information per group needed for filtering
  i = match(dataset$peptides$key_peptide_group, dataset$dt_pep_group$key_peptide_group)
  # for QC figures, consider in each group those peptides that have 1 identification (or 2 if DIA) and 2+ quantifications

  # filter
  dataset$peptides[ , colname] = dataset$peptides$intensity
  if(disregard_exclude_samples) {
    # update: if criterium is larger than group-size, limit to group-size (eg; nquant=3 while only 2 samples in group)
    dataset$peptides[!is.finite(i) | dataset$dt_pep_group$nquant_noexclude[i] < pmin(dataset$dt_pep_group$size_noexclude[i], nquant) | dataset$dt_pep_group$ndetect_noexclude[i] < pmin(dataset$dt_pep_group$size_noexclude[i], ndetect), colname] <- NA
    # dataset$peptides[is.na(i) | dataset$dt_pep_group$nquant_noexclude[i] < nquant | dataset$dt_pep_group$ndetect_noexclude[i] < ndetect, colname] <- NA
  } else {
    # analogous to above, but filter on all samples
    dataset$peptides[!is.finite(i) | dataset$dt_pep_group$nquant[i] < pmin(dataset$dt_pep_group$size[i], nquant) | dataset$dt_pep_group$ndetect[i] < pmin(dataset$dt_pep_group$size[i], ndetect), colname] <- NA
    # dataset$peptides[is.na(i) | dataset$dt_pep_group$nquant[i] < nquant | dataset$dt_pep_group$ndetect[i] < ndetect, colname] <- NA
  }

  # normalize
  dataset = normalize_peptide_intensity_column(dataset, col_intensity = colname, norm_algorithm = norm_algorithm, rollup_algorithm = rollup_algorithm)
  # dataset$peptides[ , colname] = normalize_intensities(data.table::data.table(dataset$peptides %>% select(key_peptide_sample, key_peptide, key_protein, key_sample, key_group, intensity = !!colname)), norm_algorithm = norm_algorithm)
  return(dataset)
}



#' Apply an array of normalization algorithms to some intensity column in the peptide table
#'
#' @param dataset your dataset
#' @param col_intensity intensity column in dataset$peptides tibble
#' @param norm_algorithm array of normalization algorithm character IDs
#' @param rollup_algorithm the algorithm for combining peptides to proteins as used in normalization algorithms that require a priori rollup from peptides to a protein-level abundance matrix (e.g. modebetween_protein). Refer to \code{\link{rollup_pep2prot}} function documentation for available options and a brief description of each
#' @importFrom data.table data.table is.data.table setorder dcast
#' @export
normalize_peptide_intensity_column = function(dataset, col_intensity, norm_algorithm, rollup_algorithm ) {
  if(length(col_intensity) != 1 || !is.character(col_intensity) || ! col_intensity %in% colnames(dataset$peptides)) {
    append_log(sprintf("invalid dataset$peptides column provided to normalize_peptide_intensity_column(). Provided value: %s", paste(col_intensity, collapse = ",")), type = "error")
  }
  if(!all(is.character(norm_algorithm)) || !any(norm_algorithm != "")) {
    append_log(sprintf("invalid normalization algorithm provided to normalize_peptide_intensity_column(): must be (and array of) character type and contain at least one non-empty. Provided value: %s", paste(norm_algorithm, collapse = ",")), type = "error")
  }

  # remove NA values to ensure we remove empty rows and columns
  tib_subset = dataset$peptides %>% filter(!is.na(!!as.symbol(col_intensity)))

  # only check for total absence of values, for now
  if(nrow(tib_subset) == 0) {
    return(dataset)
  }

  # convert to wide. this is the matrix we need to normalize
  tibw_subset = tib_subset %>%
    pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = !!col_intensity, values_fill = NA)

  # tibble to numeric matrix
  mat = as.matrix(tibw_subset %>% select(-peptide_id))

  # group names
  m_groups = dataset$samples$group[match(colnames(mat), dataset$samples$sample_id)] # character array, enforced upstream by sample_metadata_sort_and_filter()
  # peptide-to-protein grouping
  m_protein = tib_subset$protein_id[match(tibw_subset$peptide_id, tib_subset$peptide_id)]

  # actual normalization; apply all algorithms iteratively
  for(alg in norm_algorithm) {
    mat = normalize_matrix(x_as_log2 = mat, algorithm = alg, group_by_cols = m_groups, group_by_rows = m_protein, rollup_algorithm = rollup_algorithm)
  }

  # back to long-format
  tib_result = as_tibble(mat) %>%
    # matrix to tibble, and add peptide ID back
    add_column(peptide_id = tibw_subset$peptide_id) %>%
    # to long format
    pivot_longer(cols = -"peptide_id", names_to = "sample_id", values_to = "intensity_temp") %>%
    # remove NA from normalized data matrix
    filter(is.finite(intensity_temp)) %>%
    rename({{col_intensity}} := intensity_temp)


  # importantly, now that we allow normalization functions to impute, this function can actually EXPAND the peptide tibble
  # e.g. peptide P in sample S was never observed in the input data (not detect, not MBR) -> impute -> must expand table


  # FULL JOIN original peptide tables and our normalization results from this function
  dataset$peptides = dataset$peptides %>%
    select(-one_of(!!col_intensity)) %>%
    full_join(tib_result, by = c("peptide_id", "sample_id"))

  # carry over other columns from peptide_id. e.g. protein_id, sequence_plain, etc. set decoy and detect columns to boolean
  dataset = dataset_transfer_peptide_properties(dataset)

  # if any imputation was performed, invalidate cache
  if(anyNA(dataset$peptides$intensity)) {
    dataset = invalidate_cache(dataset)
  }

  return(dataset)
  # debug, show imputed rows; tmp %>% filter(is.na(intensity)) %>% select(!starts_with("key_"))
}



#' carries over protein_id, sequence_plain and sequence_modified in rows where protein_id is NA
#'
#' all NA values in detect and isdecoy columns will be set to FALSE
#'
#' @param dataset input dataset
#' @importFrom data.table chmatch
dataset_transfer_peptide_properties = function(dataset) {
  stopifnot(!is.na(dataset$peptides$peptide_id) & !is.na(dataset$peptides$sample_id))
  # rows that need updating; protein_id is NA. We here assume that rows with a protein_id have sequence data
  rows_todo = is.na(dataset$peptides$protein_id)
  if(!any(rows_todo)) {
    return(dataset)
  }

  # reference table
  tib_ref = dataset$peptides %>% filter(!rows_todo)
  # index from entire table to reference. use data.table::chmatch as it's a bit faster than regular match
  i = data.table::chmatch(dataset$peptides$peptide_id, tib_ref$peptide_id)

  # just blanked overwrite, faster than subsetting all the time
  dataset$peptides$protein_id = tib_ref$protein_id[i]
  dataset$peptides$sequence_plain = tib_ref$sequence_plain[i]
  dataset$peptides$sequence_modified = tib_ref$sequence_modified[i]

  # optional columns to set to a default value
  if("detect" %in% colnames(dataset$peptides)) {
    dataset$peptides$detect[is.na(dataset$peptides$detect)] = FALSE
  }
  if("isdecoy" %in% colnames(dataset$peptides)) {
    dataset$peptides$isdecoy[is.na(dataset$peptides$isdecoy)] = FALSE
  }
  return(dataset)
}
