
#' wrapper function for peptide-to-protein rollup
#'
#' @param tib peptide tibble in long-format. missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating log2 scaling
#' @param algo_rollup type of rollup to perform. options: maxlfq or maxlfq_iq for MaxLFQ rollup using iq package, maxlfq_diann for MaxLFQ rollup using diann package. Any other value performs rollup by peptide summation.
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @export
rollup_pep2prot = function(tib, intensity_is_log2 = TRUE, algo_rollup, return_as_matrix = FALSE) {
  if(algo_rollup == "maxlfq" || algo_rollup == "maxlfq_iq") {
    return(rollup_pep2prot_maxlfq(tib, intensity_is_log2 = intensity_is_log2, implementation = "iq", return_as_matrix = return_as_matrix))
  }
  if(algo_rollup == "maxlfq_diann") {
    return(rollup_pep2prot_maxlfq(tib, intensity_is_log2 = intensity_is_log2, implementation = "diann", return_as_matrix = return_as_matrix))
  }
  return(rollup_pep2prot_summation(tib, intensity_is_log2 = intensity_is_log2, return_as_matrix = return_as_matrix))
}



#' rollup by plain peptide intensity summation
#'
#' @param tib peptide tibble in long-format. missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating log2 scaling
#' @param return_as_matrix instead of long-format tibble, return a protein*sample matrix (missing values set to NA)
#' @importFrom tidyr pivot_wider pivot_longer
#' @export
rollup_pep2prot_summation = function(tib, intensity_is_log2 = TRUE, return_as_matrix = FALSE) {
  stopifnot(c("peptide_id", "protein_id", "sample_id", "intensity") %in% colnames(tib))
  stopifnot(intensity_is_log2 %in% c(T,F))
  stopifnot(return_as_matrix %in% c(T,F))
  start_time = Sys.time()

  ### format intensities
  tib = tib %>% filter(is.finite(intensity))
  # undo log2 transformation before summation
  if(intensity_is_log2) {
    tib$intensity = 2^tib$intensity
  }

  # convert to matrix and aggregate
  tibw = tib %>% filter(is.finite(intensity)) %>% tidyr::pivot_wider(id_cols = c("peptide_id", "protein_id"), names_from = "sample_id", values_from = "intensity", values_fill = list(intensity=0))
  df_prot = aggregate.data.frame(tibw %>% select(-peptide_id, -protein_id), by = list(protein_id = tibw$protein_id), FUN = sum, na.rm=T)


  # if a matrix was requested we're done
  if(return_as_matrix) {
    mat_prot = as_matrix_except_first_column(df_prot)
    mat_prot[!is.finite(mat_prot) | mat_prot <= 0] = NA
    # sync log transformation with input data (currently log2 transformed)
    if(intensity_is_log2) {
      mat_prot = log2(mat_prot)
    }
    return(mat_prot)
  }

  # back to long-format & remove zeros
  tib_result = as_tibble(df_prot) %>% tidyr::pivot_longer(cols = -protein_id, names_to = "sample_id", values_to = "intensity") %>%
    filter(is.finite(intensity) & intensity > 0)

  # sync log transformation with input data (currently not log transformed)
  if(intensity_is_log2) {
    tib_result = tib_result %>% mutate(intensity = log2(intensity))
  }

  append_log_timestamp("peptide to protein rollup by summation of intensities", start_time) # maxlfq is pretty slow, so print timer to console so the users are aware this is a time consuming step
  return(tib_result)
  # one-liner for reference (but slower); expr_prot = as_tibble(2^x_as_log2) %>% replace(is.na(.), 0) %>% add_column(protein_id=grp_var) %>% group_by(protein_id) %>% summarise_all(sum) %>% replace(.==0, NA)
}



#' rollup by MaxLFQ algorithm
#'
#' Available implementations:
#' iq R package: https://github.com/tvpham/iq/
#' DIA-NN R package: https://github.com/vdemichev/DiaNN/
#'
#' @param tib peptide tibble in long-format. missing values should be set NA, not zero. required columns: peptide_id, protein_id, sample_id, intensity
#' @param intensity_is_log2 boolean indicating log2 scaling
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
  # diann maxlfq implementation will crash on zero values
  # we simply threshold at 1, _assuming_ upstream data providers (eg; intensity values from user input dataset) are well above zero and any incidental value < 1 is caused by normalization of values that were already near zero
  tib$intensity[is.finite(tib$intensity) & tib$intensity < 0] = 0

  if(implementation == "diann") {
    # use the MaxLFQ implementation provided by the DIA-NN team. output is a plain matrix
    tib$intensity[is.finite(tib$intensity) & tib$intensity < 0.1] = 0.1 # tested diann maxlfq implementation crashed on some datasets containing zero values, so the thresholded values set to zero upstream are here moved to 0.1
    mat_prot = diann::diann_maxlfq(tib, sample.header = "sample_id", group.header="protein_id", id.header = "peptide_id", quantity.header = "intensity")
  }
  if(implementation == "iq") {
    # use the MaxLFQ implementation provided by the iq R package. output contains a plain matrix at $estimate
    invisible(capture.output( mat_prot <- iq::fast_MaxLFQ( list(protein_list = tib$protein_id, sample_list = tib$sample_id, id = tib$peptide_id, quant = tib$intensity) )$estimate  )) # note, must use <- notation within capture.output
  }

  # sync log transformation with input data (currently log2 transformed)
  if( ! intensity_is_log2) {
    mat_prot = 2^mat_prot
  }

  # if a matrix was requested we're done
  if(return_as_matrix) {
    mat_prot[!is.finite(mat_prot) | mat_prot <= 0] = NA
    result = mat_prot
  } else {
    # back to long-format & remove zeros
    result = as_tibble(mat_prot) %>% add_column(protein_id = rownames(mat_prot)) %>% tidyr::pivot_longer(cols = -protein_id, names_to = "sample_id", values_to = "intensity") %>%
      filter(is.finite(intensity) & intensity > 0)
  }

  append_log_timestamp(sprintf("peptide to protein rollup with MaxLFQ (implementation: %s)", implementation), start_time) # maxlfq is pretty slow, so print timer to console so the users are aware this is a time consuming step
  return(result)

  # debug;
  # tib_prot_summ = rollup_pep2prot_summation(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T)
  # tib_prot1 = rollup_pep2prot_maxlfq(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T, implementation = "iq")
  # tib_prot2 = rollup_pep2prot_maxlfq(tib=dataset$peptides %>% select(peptide_id, protein_id, sample_id, intensity = intensity_all_group), intensity_is_log2 = T, implementation = "diann")
  # tib_plot = tib_prot1 %>% rename(intensity_iq=intensity) %>% left_join(tib_prot2 %>% rename(intensity_diann=intensity))
  # plot(tib_plot$intensity_iq, tib_plot$intensity_diann, main=sd(tib_plot$intensity_iq - tib_plot$intensity_diann))
  # table(tib_plot$intensity_iq == tib_plot$intensity_diann)
  #
  # tib_plot = tib_prot1 %>% rename(intensity_iq=intensity) %>% left_join(tib_prot_summ %>% rename(intensity_summ=intensity))
  # plot(tib_plot$intensity_iq, tib_plot$intensity_summ, main=sd(tib_plot$intensity_iq - tib_plot$intensity_summ))
  #
  # tib %>% filter(protein_id == "P55011") %>% filter(is.finite(intensity))
  # tib_prot_summ %>% filter(protein_id == "P55011")
  # tib_prot1 %>% filter(protein_id == "P55011")
  # tib_prot2 %>% filter(protein_id == "P55011")
}


