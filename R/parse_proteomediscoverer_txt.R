####
#### ProteomeDiscoverer bug;
####
#
## PSMs.tsv (looks fine)
# Annotated Sequence	Modifications
# [K].AAGLATmISTmRPDIDNmDEYVR.[N]	M7(Oxidation); M11(Oxidation); M18(Oxidation)
# [K].AAGLATMISTmRPDIDNmDEYVR.[N]	M11(Oxidation); M18(Oxidation)
# [K].AAGLATmISTmRPDIDNMDEYVR.[N]	M7(Oxidation); M11(Oxidation)
# [K].AAGLATMISTmRPDIDNMDEYVR.[N]	M11(Oxidation)
# [K].AAGLATMISTmRPDIDNMDEYVR.[N]	M11(Oxidation)
#
## PeptideGroups.tsv (malformed modification annotations, see second row)
# Annotated Sequence	Modifications
# [K].AAGLATMISTMRPDIDNMDEYVR.[N]	1xOxidation [M11]
# [K].AAGLATMISTMRPDIDNMDEYVR.[N]	2xOxidation [M11; M]
# [K].AAGLATMISTMRPDIDNMDEYVR.[N]	3xOxidation [M7; M11; M18]
#
#
# 1xAcetyl [N-Term]; 1xOxidation [M2]



#' Import a label-free DDA proteomics dataset from ProteomeDiscoverer, experimental feature !
#'
#' Input data must contain Percolator PEP scores for each PSM, so after search engine (eg; Sequest HT) make sure to connect the Percolator node.
#'
#' Example PD workflow:
#' Processing Step: PWF_QE_Precursor_Quan_and_LFQ_SequestHT_Percolator
#' Consensus Step: CWF_Comprehensive_Enhanced Annotation_LFQ_and_Precursor_Quan
#' Consensus Step: add the "result exporter" (drag&drop from side panel to bottom panel)
#'
#' Optionally, you can relax the output filter criteria somewhat, thereby relying on MS-DAP filtering a bit more and potentially gaining some extra (MBR)hits, as follows;
#' Consensus step -->> "peptide and protein filter" -->> Peptide Confidence At Least -->> change to medium
#' While you're making changes there, you can problably set "Remove Peptides Without Proteins" to "True"
#'
#'
#'
#' @param filename full path to the _PSMs.txt file
#'
#' @importFrom data.table fread chmatch
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest
#' @importFrom stringr str_sub
#' @export
#'
import_dataset_proteomediscoverer_txt = function(filename) {
  # debug;
  # dataset = import_dataset_proteomediscoverer_txt(filename = "C:/DATA/PXD007683/pd/PXD007683/pdout_PSMs.txt")
  # x = dataset$peptides %>% pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = "intensity")
  # plot(x$a05191, x$a05192, main = sprintf("R^2=%.4f", cor(x$a05191, x$a05192, use = "pairwise.complete.obs", method = "pearson")^2))

  reset_log()
  append_log("ProteomeDiscoverer data import is a new feature, please consult the documentations for limitations / work-in-progress", type = "warning")

  # first, check if input file exists
  check_parameter_is_string(filename)
  if (!file.exists(filename)) {
    append_log(paste("input file does not exist:", filename), type = "error")
  }


  tib_psm = tibble::as_tibble(data.table::fread(filename))

  cols_required = c("Percolator PEP", "Confidence", "# Protein Groups", "Annotated Sequence", "Modifications", "Charge", "Spectrum File", "RT [min]", "Precursor Abundance", "Master Protein Accessions")
  cols_missing = setdiff(cols_required, colnames(tib_psm))
  if(length(cols_missing) > 0) {
    append_log(paste("ProteomeDiscoverer input file must be *_PSMs.txt from a label-free quantification workflow that includes Percolator. The following columns are missing from provided input file:", paste(cols_missing, collapse = ", ")), type = "error")
  }

  tib_psm$`Percolator PEP` = suppressWarnings(as.numeric(tib_psm$`Percolator PEP`))
  tib_psm$`# Protein Groups` = suppressWarnings(as.numeric(tib_psm$`# Protein Groups`))
  tib_psm$`RT [min]` = suppressWarnings(as.numeric(tib_psm$`RT [min]`))
  tib_psm$`Precursor Abundance` = suppressWarnings(as.numeric(tib_psm$`Precursor Abundance`))

  # throw away everything not even close to significant
  tib_psm = tib_psm %>%
    filter(is.finite(`Percolator PEP`) & `Percolator PEP` <= 0.1 &
             Confidence %in% c("High", "Medium") &
             is.finite(`RT [min]`) &
             `Master Protein Accessions` != "" &
             `Spectrum File` != "" &
             `Annotated Sequence` != "" &
             is.finite(`Precursor Abundance`) &
             # throw away all PSM/peptides that map to multiple proteinsgroups
            `# Protein Groups` == 1) %>%
    arrange(`Percolator PEP`)

  ## create some unique keys
  tib_psm$peptide_id = paste0(tib_psm$`Annotated Sequence`, "_", tib_psm$Modifications, "_", tib_psm$Charge)
  tib_psm$peptide_x_sample = paste0(tib_psm$peptide_id, "***", tib_psm$`Spectrum File`)
  # debug; head(tib_psm$peptide_id); head(tib_psm$peptide_x_sample)

  ## if there are multiple PSM entries for a precursor in some sample, we keep the most significant (best MS/MS score)
  tib_psm = tib_psm %>% distinct(peptide_x_sample, .keep_all = TRUE)

  ## suppose that we removed a redundant PSM which was significant AND more abundant, prefer that
  # simplest implementation is to do a reverse lookup and overwrite abundances.
  # to re-iterate; for _all_ significant peptides we prefer the most abundant and for those not significant we prefer the most confident
  tib_psm_signif = tib_psm %>%
    filter(`Percolator PEP` <= 0.01 & Confidence %in% c("High", "Medium")) %>%
    arrange(desc(`Precursor Abundance`))

  i = data.table::chmatch(tib_psm_signif$peptide_x_sample, tib_psm$peptide_x_sample)
  tib_psm$`Precursor Abundance`[i] = tib_psm_signif$`Precursor Abundance`
  tib_psm$`RT [min]`[i] = tib_psm_signif$`RT [min]`
  rm(tib_psm_signif)

  # hist(log10(tib_psm$`Percolator PEP`)); hist(log10(tib_psm$`Precursor Abundance`))

  ## subset and format tibble for compatability with out pipeline
  tib_result = tib_psm %>%
    select(sample_id = `Spectrum File`, sequence_input = `Annotated Sequence`, peptide_id, protein_id=`Master Protein Accessions`, charge = Charge, Modifications, intensity = `Precursor Abundance`, rt = `RT [min]`, confidence = `Percolator PEP`) %>%
    filter(is.finite(rt) & is.finite(intensity) & intensity > 0)

  tryCatch({
    tib_result$sequence_key = paste0(tib_result$sequence_input, "_", tib_result$Modifications)

    # cleanup raw file names
    # regex from msdap::import_dataset_in_long_format()
    tib_result$sample_id = gsub("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia)(\\.gz|\\.dia){0,1}$)", "", tib_result$sample_id, ignore.case=T)


    ### cleanup sequences
    # subset of all unique sequences (so "expensive" string operations are limited)
    tib_seq = tib_result %>% select(key = sequence_key, input = sequence_input, modifications = Modifications) %>% distinct(key, .keep_all = TRUE)

    # plain sequence = string everything outside the brackets and dots + uppercase
    tib_seq$sequence_plain = tib_seq$sequence_modified = toupper(gsub("(^\\.)|(^(.*\\]\\.))|(\\.$)|(\\.\\[.*)", "", tib_seq$input))

    ## pretty-print modified sequences
    ## example; "N-Term(Prot)(Met-loss+Acetyl); M16(Oxidation); C17(Carbamidomethyl)"
    # subset only sequences with a modification for efficiency
    # tib_seq = tib_seq %>% filter(modifications != "")
    tib_seq$index_before_unnest = 1:nrow(tib_seq)
    # split multiple mods
    tib_seq$modlist = strsplit(tib_seq$modifications, ";", fixed = T)


    tib_seq_unnested = tib_seq %>% tidyr::unnest(cols = modlist)
    # extract name and position
    tib_seq_unnested$name = gsub(".*(\\([^(]+\\))[^(]*", "\\1", tib_seq_unnested$modlist)
    tib_seq_unnested$pos = gsub("^\\D*(\\d+)\\D.*", "\\1", tib_seq_unnested$modlist)
    tib_seq_unnested$pos[grepl("N-Term", tib_seq_unnested$modlist, fixed = T)] = "0"
    tib_seq_unnested$pos = suppressWarnings(as.integer(tib_seq_unnested$pos))
    rows_fail = !is.finite(tib_seq_unnested$pos)
    if(any(rows_fail)) {
      tib_seq_unnested$pos[rows_fail] = nchar(tib_seq_unnested$sequence_plain[rows_fail])
    }

    # sort
    tib_seq_unnested = tib_seq_unnested %>% arrange(desc(pos))

    # pretty-printed sequences
    for(j in 1:nrow(tib_seq_unnested)) { #j=5
      i = tib_seq_unnested$index_before_unnest[j]
      # tib_seq_unnested[j,]
      tib_seq$sequence_modified[i] = paste0(stringr::str_sub(tib_seq$sequence_modified[i], 0, tib_seq_unnested$pos[j]), tib_seq_unnested$name[j], stringr::str_sub(tib_seq$sequence_modified[i], tib_seq_unnested$pos[j] + 1))
    }

    # map to result tibble
    i = data.table::chmatch(tib_result$sequence_key, tib_seq$key)
    tib_result$sequence_plain = tib_seq$sequence_plain[i]
    tib_result$sequence_modified = tib_seq$sequence_modified[i]

    # use the pretty-printed modified sequences to re-define peptide_id, then cleanup all temporary columns
    tib_result$peptide_id = paste0(tib_result$sequence_modified, "_", tib_result$charge)
  }, error = function(e) {
    tib_result = tib_psm %>%
      select(sample_id = `Spectrum File`, sequence_plain = `Annotated Sequence`, peptide_id, protein_id=`Master Protein Accessions`, charge = Charge, Modifications, intensity = `Precursor Abundance`, rt = `RT [min]`, confidence = `Percolator PEP`) %>%
      filter(is.finite(rt) & is.finite(intensity) & intensity > 0) %>%
      mutate(sequence_modified = paste0(sequence_plain,  "_", Modifications))
  })

  tib_result$sequence_input = NULL
  tib_result$sequence_key = NULL
  tib_result$Modifications = NULL
  # debug; print(tib_result, n=99)

  # ensure downstream compatability
  tib_result$isdecoy = F
  tib_result$intensity = log2(tib_result$intensity)
  tib_result$intensity[!is.na(tib_result$intensity) & tib_result$intensity < 0] = 0 # note; we already removed zero intensity values when importing. here, we threshold extremely low values
  tib_result$detect = tib_result$confidence <= 0.01 # peptide-level FDR cutoff


  tib_result = peptides_collapse_by_sequence(tib_result, prop_peptide = "sequence_modified")

  log_peptide_tibble_pep_prot_counts(tib_result)
  return(list(peptides=tibble_peptides_reorder(tib_result), proteins=empty_protein_tibble(tib_result), acquisition_mode = "dda"))
}







# df = data.frame(x=(tib_pepgroup$`Abundances (Scaled): F1: Sample`), y=(tib_pepgroup$`Abundances (Scaled): F2: Sample`), stringsAsFactors = F)
# df = df[is.finite(df$x) & is.finite(df$y), ]
# plot(df)
# hist(df$x)
# hist(df$y)
#
# plot(df$x, df$y)
# hist(df$x - df$y)
# head(df)

#
# file_pepgroup = "C:/DATA/PXD007683/pd/PXD007683/testfile_PeptideGroups.txt"
# tib_pepgroup = tibble::as_tibble(data.table::fread(file_pepgroup))
#
# # only peptides that map to 1 peptideGroup
# tib_pepgroup = tib_pepgroup %>% filter(`# Protein Groups` == 1)
#
#
# ### extract peptide intensities
# col_intensity = grep("^Abundances \\(Scaled\\): ", colnames(tib_pepgroup), value=T, ignore.case = T)
#
# # subset data + parse sample names
# tib_temp = tib_pepgroup %>% select(all_of(c("Annotated Sequence", col_intensity)))
# colnames(tib_temp) = gsub(".*: *(F\\d+) *:.*", "\\1", colnames(tib_temp), ignore.case = T)
#
# # wide-to-long
# tib = tidyr::pivot_longer(tib_temp, cols = -c("Annotated Sequence"), names_to = "sample_id", values_to = "intensity") %>%
#   mutate(intensity = suppressWarnings(as.numeric(intensity))) %>%
#   filter(is.finite(intensity))
#
# # add peptide info
# tib = tib %>% left_join(tib_pepgroup %>% select(protein_id = `Master Protein Accessions`, `Annotated Sequence`, modifications = Modifications), by="Annotated Sequence")
#
#
# ## modified sequences
# tib_seq = tib_pepgroup %>% select(`Annotated Sequence`, Modifications)
# tib_seq$sequence_modified = tib_seq$sequence_plain = gsub("(^\\.)|(^(.*\\]\\.))|(\\.$)|(\\.\\[.*)", "", tib_seq$`Annotated Sequence`)
#
#
#
#
# ## inject modifications by manual string manipulation
# decorate_modified_sequences = function(tib_seq) {
#   ##### while below code may seem unlogical / not-so-clean, we choose relatively fast functions over cleanly formatting in intermediate steps + maximize vectorization of all string manipulation
#   # regex is life
#
#   # example: 1xCarbamidomethyl [C18]; 2xOxidation [M11; M19]
#   i = which(!is.na(tib_seq$Modifications) & !tib_seq$Modifications %in% c("", "null", "na"))
#   # find split points (multiple distinct modifications)
#   mods = gsub("\\] *; *(\\d)", "]@@@@\\1", tib_seq$Modifications[i])
#   mods = strsplit(mods, "@@@@", fixed = T)
#   x = tibble(index = rep(i, lengths(mods)), mods = unlist(mods, use.names = F))
#   # x[x$index == i[997],] # debug; check an example
#
#   # split entries again (eg; 2xOxidation [M11; M19])
#   x$name = gsub("( *\\d+ *x *)| *\\[.*", "", x$mods)
#   x$pos = strsplit(x$mods, ";", fixed = T)
#   x = x %>% tidyr::unnest(cols = pos)
#
#   # cleanup the position
#   x$pos = sub("(^|.*\\D)(\\d+)\\D*$", "\\2", x$pos)
#
#   ### so now, x is a tibble that holds;
#   # index = row in the tib_seq table
#   # name = pretty-print name of the modification
#   # pos = character position in tib_seq$sequence_plain (at given index) where the modification should be inserted
#
#   # sort the table such that we start inserting at the end; that way we don't invalidate string positions (nor do we need to keep track of padding)
#   # x = x %>% mutate(pos = as.integer(pos)) %>% arrange(desc(pos))
#   x[!is.finite(x$pos),]
#   tib_seq[120, ]
#
#   for(j in 1:nrow(x)) { #j=1
#     tib_seq$Modifications[x$index[j]] = paste0(stringr::str_sub(tib_seq$Modifications[x$index[j]], 0, x$pos[j]), "(", x$name[j], ")", stringr::str_sub(tib_seq$Modifications[x$index[j]], x$pos[j] + 1))
#   }
#
#   return(result)
# }
#
#
#
# # Annotated Sequence
# # Modifications
# # Master Protein Accessions
# # Abundances (Scaled): F1: Sample	Abundances (Scaled): F2: Sample
# # Qvality PEP # only for first-pass filter
#
#
# # example;
# # [K].FLNNTYENPFMNASGVHCMTTQELDELANSK.[A]	1xCarbamidomethyl [C18]; 2xOxidation [M11; M19]
# # P00924; P00925
#
#
#
# #### PSM
# #
# # RT [min]
# # Spectrum File
# # Identifying Node
# # Annotated Sequence
# # PSM Ambiguity
# # `# Protein Groups`
# # Master Protein Accessions
# # Precursor Abundance
# # Charge
# # m/z [Da]
#
#
# # select(peptide_id, protein_id, sample_id = raw_file, sequence_plain, sequence_modified, charge, mz, intensity, confidence = pep, isdecoy, rt = calibrated_rt, detect)
#
