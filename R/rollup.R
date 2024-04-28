
#' wrapper function for peptide-to-protein rollup
#'
#' converts a table with peptide-level abundance values to protein-level data by combining data from respective peptides per protein (per sample).
#'
#' @param tib peptide tibble in long-format. missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating if input peptide intensities are log2 scaled
#' @param rollup_algorithm algorithm for combining peptides to proteins. Options: 'maxlfq', 'sum', 'tukey_median'. The former applies the MaxLFQ algorithm, the second employs the 'classic' strategy of summing all peptides per protein (non-log values), the latter uses Tukey's median polish algorithm. Alternatively, one can specify "maxlfq_diann" to use the DIA-NN R package's implementation for MaxLFQ instead of the faster implementation from the iq R package (default when using 'maxlfq')
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @export
rollup_pep2prot = function(tib, intensity_is_log2 = TRUE, rollup_algorithm, return_as_matrix = FALSE) {
  if (length(rollup_algorithm) != 1 || ! rollup_algorithm %in% c("sum", "maxlfq_diann", "maxlfq", "tukey_median")) {
    append_log("rollup_algorithm parameter only supports 'sum', 'maxlfq', 'maxlfq_diann' and 'tukey_median'", type = "error")
  }

  if(rollup_algorithm == "maxlfq") {
    return(rollup_pep2prot_maxlfq(tib, intensity_is_log2 = intensity_is_log2, implementation = "iq", return_as_matrix = return_as_matrix))
  }
  if(rollup_algorithm == "maxlfq_diann") {
    return(rollup_pep2prot_maxlfq(tib, intensity_is_log2 = intensity_is_log2, implementation = "diann", return_as_matrix = return_as_matrix))
  }
  if(rollup_algorithm == "sum") {
    return(rollup_pep2prot_summation(tib, intensity_is_log2 = intensity_is_log2, return_as_matrix = return_as_matrix))
  }
  if(rollup_algorithm == "tukey_median") {
    return(rollup_pep2prot_tmp(tib, intensity_is_log2 = intensity_is_log2, return_as_matrix = return_as_matrix))
  }

  append_log(paste("unknown rollup_algorithm parameter;", rollup_algorithm), type = "error")
}



#' rollup a peptide tibble by summation of all intensity values per protein
#'
#' @param tib peptide tibble in long-format. Required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating if input peptide intensities are log2 scaled
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @importFrom tidyr pivot_wider pivot_longer
#' @export
rollup_pep2prot_summation = function(tib, intensity_is_log2 = TRUE, return_as_matrix = FALSE) {
  stopifnot(is.data.frame(tib) && length(tib) > 0)
  stopifnot(c("peptide_id", "protein_id", "sample_id", "intensity") %in% colnames(tib))
  stopifnot(intensity_is_log2 %in% c(T,F))
  stopifnot(return_as_matrix %in% c(T,F))
  start_time = Sys.time()

  # sum using non-log intensities
  if(intensity_is_log2) {
    tib$intensity = 2^tib$intensity
  }

  # convert peptide tibble to wide-format for efficient aggregation downstream. Remove non-finite intensities, then fill missing values with 0
  tibw = tib %>%
    filter(is.finite(intensity)) %>%
    tidyr::pivot_wider(id_cols = c("peptide_id", "protein_id"), names_from = "sample_id", values_from = "intensity", values_fill = list(intensity=0))

  # rollup
  df_prot = stats::aggregate.data.frame(tibw %>% select(-peptide_id, -protein_id), by = list(protein_id = tibw$protein_id), FUN = sum, na.rm=T)

  if(return_as_matrix) {
    # if output should be a matrix...
    result = as_matrix_except_first_column(df_prot)
    result[!is.finite(result) | result <= 0] = NA
    # protein intensity values are not log2 atm. If input was, sync by applying log2 transformation
    if(intensity_is_log2) {
      result = log2(result)
    }
  } else {
    # otherwise, convert from wide-format to long-format
    result = as_tibble(df_prot) %>%
      tidyr::pivot_longer(cols = -protein_id, names_to = "sample_id", values_to = "intensity") %>%
      filter(is.finite(intensity) & intensity > 0)
    # protein intensity values are not log2 atm. If input was, sync by applying log2 transformation
    if(intensity_is_log2) {
      result$intensity = log2(result$intensity)
    }
  }

  append_log_timestamp("peptide to protein rollup by summation", start_time)
  return(result)
  # one-liner for reference (but slower); expr_prot = as_tibble(2^x_as_log2) %>% replace(is.na(.), 0) %>% add_column(protein_id=grp_var) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA)
}



#' rollup a peptide tibble by applying Tukey's median polish algorithm to all peptides per protein
#'
#' @param tib peptide tibble in long-format. Missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating if input peptide intensities are log2 scaled
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @importFrom tidyr pivot_wider pivot_longer
#' @export
rollup_pep2prot_tmp = function(tib, intensity_is_log2 = TRUE, return_as_matrix = FALSE) {
  stopifnot(c("peptide_id", "protein_id", "sample_id", "intensity") %in% colnames(tib))
  stopifnot(intensity_is_log2 %in% c(T,F))
  stopifnot(return_as_matrix %in% c(T,F))
  start_time = Sys.time()

  # for TMP, use log intensities
  if(!intensity_is_log2) {
    tib$intensity = suppressWarnings(log2(tib$intensity)) # non-finites will be filtered in next step
  }

  # peptide*sample matrix of log2 intensity values. Missing values are NA
  tibw = tib %>%
    filter(is.finite(intensity)) %>%
    tidyr::pivot_wider(id_cols = c("peptide_id", "protein_id"), names_from = "sample_id", values_from = "intensity", values_fill = list(intensity=NA))
  tibw_mat = as.matrix(tibw %>% select(-peptide_id, -protein_id))

  # lookup table from protein_id to index in peptide wide-format table `tibw`
  tibw_lookup = tibw %>% add_column(index = 1:nrow(tibw)) %>% group_by(protein_id) %>% summarise(index = list(index), .groups = "drop")

  # iterate unique protein IDs
  mat_prot = matrix(0, nrow = nrow(tibw_lookup), ncol = ncol(tibw_mat), dimnames = list(tibw_lookup$protein_id, colnames(tibw_mat)) )
  for(i in 1:nrow(mat_prot)) { #i=1
    rows = tibw_lookup$index[[i]] # rows in peptide matrix `tibw` for current protein i
    if(length(rows) == 1) {
      # only 1 peptide, nothing to do
      mat_prot[i,] = tibw_mat[rows,]
    } else {
      # use base R implementation for median polish
      # to go faster, parallelise and use Rcpp implementation of TMP algorithm
      tmp = suppressWarnings(stats::medpolish(tibw_mat[rows,], na.rm = TRUE, trace.iter = FALSE))
      mat_prot[i,] = tmp$overall + tmp$col
    }
  }
  mat_prot[!is.finite(mat_prot) | mat_prot <= 0] = NA

  # `mat_prot` is log2 prior to this line. If input wasn't, sync by reversing log2 transformation
  if(!intensity_is_log2) {
    mat_prot = 2^mat_prot
  }

  result = mat_prot
  if(!return_as_matrix) {
    # from matrix to long-format tibble
    result = as_tibble(mat_prot) %>%
      add_column(protein_id = rownames(mat_prot), .before = 1) %>%
      tidyr::pivot_longer(cols = -protein_id, names_to = "sample_id", values_to = "intensity") %>%
      filter(is.finite(intensity))
  }

  append_log_timestamp("peptide to protein rollup by Tukey's median polish", start_time)
  return(result)
}



#' rollup by MaxLFQ algorithm
#'
#' Available implementations:
#' iq R package: https://github.com/tvpham/iq/
#' DIA-NN R package: https://github.com/vdemichev/DiaNN/
#'
#' @param tib peptide tibble in long-format. missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating if input peptide intensities are log2 scaled
#' @param implementation implementation to use. available options: "iq", "diann"
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @importFrom diann diann_maxlfq
#' @importFrom iq fast_MaxLFQ
#' @export
rollup_pep2prot_maxlfq = function(tib, intensity_is_log2, implementation = "iq", return_as_matrix = FALSE) {
  stopifnot(c("peptide_id", "protein_id", "sample_id", "intensity") %in% colnames(tib))
  stopifnot(intensity_is_log2 %in% c(T,F))
  stopifnot(return_as_matrix %in% c(T,F))
  stopifnot(implementation %in% c("iq", "diann"))
  start_time = Sys.time()

  ### format intensities
  tib = tib %>% filter(is.finite(intensity))
  # downstream function expects log2 transformed data
  if( ! intensity_is_log2) {
    tib$intensity = log2(tib$intensity)
  }

  ### MaxLFQ
  if(implementation == "diann") {
    # use the MaxLFQ implementation provided by the DIA-NN team. output is a plain matrix
    # diann maxlfq implementation will crash on zero values
    tib$intensity[is.finite(tib$intensity) & tib$intensity < 0.1] = 0.1 # tested diann maxlfq implementation crashed on some datasets containing zero values, so the thresholded values set to zero upstream are here moved to 0.1
    mat_prot = diann::diann_maxlfq(tib, sample.header = "sample_id", group.header="protein_id", id.header = "peptide_id", quantity.header = "intensity")
  }
  if(implementation == "iq") {
    # use the MaxLFQ implementation provided by the iq R package. output contains a plain matrix at $estimate
    tib$intensity[is.finite(tib$intensity) & tib$intensity < 0] = 0
    invisible(capture.output( mat_prot <- iq::fast_MaxLFQ( list(protein_list = tib$protein_id, sample_list = tib$sample_id, id = tib$peptide_id, quant = tib$intensity) )$estimate  )) # note, must use <- notation within capture.output
  }

  ### format results
  # sync log transformation with input data (currently log2 transformed)
  if( ! intensity_is_log2) {
    mat_prot = 2^mat_prot
  }
  if(return_as_matrix) {
    # if output should be a matrix...
    mat_prot[!is.finite(mat_prot) | mat_prot <= 0] = NA
    result = mat_prot
  } else {
    # otherwise, convert from wide-format to long-format
    result = as_tibble(mat_prot) %>%
      add_column(protein_id = rownames(mat_prot)) %>%
      tidyr::pivot_longer(cols = -protein_id, names_to = "sample_id", values_to = "intensity") %>%
      filter(is.finite(intensity) & intensity > 0)
  }

  append_log_timestamp(sprintf("peptide to protein rollup with MaxLFQ (implementation: %s)", implementation), start_time)
  return(result)

  # debug;
  # f <- system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")
  # dataset = import_dataset_skyline(f, confidence_threshold = 0.01, return_decoys = F, acquisition_mode = "dia")
  # dataset = sample_metadata_custom(dataset, group_regex_array = c(A = "007|009|011", B = "008|010|012") )
  # dataset = filter_dataset(dataset, filter_min_detect = 3, norm_algorithm = c("vwmb", "modebetween_protein"), by_group = F, all_group = T, by_contrast = F)
  #
  # tib_prot_summ = rollup_pep2prot_summation(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T)
  # tib_prot1 = rollup_pep2prot_maxlfq(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T, implementation = "iq")
  # tib_prot2 = rollup_pep2prot_maxlfq(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T, implementation = "diann")
  #
  # tib_plot = tib_prot1 %>% rename(intensity_iq=intensity) %>% left_join(tib_prot2 %>% rename(intensity_diann=intensity)) %>% left_join(tib_prot_summ %>% rename(intensity_summ=intensity))
  # plot(tib_plot$intensity_iq, tib_plot$intensity_summ, main=sd(tib_plot$intensity_iq - tib_plot$intensity_summ))
  # plot(tib_plot$intensity_iq, tib_plot$intensity_diann, main=sd(tib_plot$intensity_iq - tib_plot$intensity_diann))
  #
  # dataset$peptides %>% filter(!is.na(intensity_all_group)) %>% count(sample_id, protein_id) %>% filter(n>1)
  # dataset$peptides %>% filter(protein_id == "1/sp|A6NL28|TPM3L_HUMAN" & is.finite(intensity_all_group)) %>% select(sample_id, peptide_id, intensity_all_group) %>% arrange(peptide_id, sample_id)
  # tib_plot %>% filter(protein_id == "1/sp|A6NL28|TPM3L_HUMAN")
}
