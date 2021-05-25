
#' Import a label-free proteomics dataset from EncyclopeDIA
#'
#' This function combines all 'ELIB' files from an EncyclopeDIA analysis into a dataset.
#'
#' Input files in the provided input directory should be the .elib files matching each sample (mzML/dia input file to EncyclopeDIA)
#'
#' @param path the directory that contains the EncyclopeDIA ELIB files
#' @param confidence_threshold confidence score threshold at which a peptide is considered 'identified' (target value must be lesser than or equals)
#' @param return_decoys logical indicating whether to return decoy peptides. Should be set to FALSE, and if enabled, make sure to manually remove the decoys from the peptides tibble before running the quickstart function!
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbGetQuery dbDisconnect
#' @importFrom data.table as.data.table
#' @export
import_dataset_encyclopedia = function(path, confidence_threshold = 0.01, return_decoys = FALSE) {
  reset_log()
  append_log("reading EncyclopeDIA elib ...", type = "info")

  check_parameter_is_string(path)
  if (!dir.exists(path)) {
    append_log(paste("input directory does not exist:", path), type = "error")
  }

  # find files with elib extension
  files = dir(path, "\\.elib$", full.names = T, ignore.case = T)
  if (length(files) < 2) {
    append_log(paste("There must be at least 2 .elib files at input directory:", path), type = "error")
  }

  peptides = pep2prot = NULL
  for(filename in files) {
    append_log(paste("parsing EncyclopeDIA elib file:", filename, "..."), type = "progress")
    # connect to database  +  catch errors
    db = try(DBI::dbConnect(RSQLite::SQLite(), filename, synchronous = NULL))
    if(length(db) == 0 || !DBI::dbIsValid(db) || tryCatch({DBI::dbExistsTable(db,"test"); FALSE}, error=function(x){TRUE})) {
      append_log(paste("input file is not an SQLite database (which the EncyclopeDIA elib file should be):", filename), type = "error")
    }

    # check if required data tables are present
    tbls = tryCatch(DBI::dbListTables(db), error=function(x){NULL})
    tbls_required = c("peptidequants", "peptidescores", "peptidetoprotein")
    if(!all(tbls_required %in% tbls)) {
      append_log(paste("EncyclopeDIA SQLite database must contain these tables:", paste(tbls_required, collapse=", ")), type = "error")
    }

    # extract data from SQLite database
    p = as_tibble(DBI::dbGetQuery(db, 'SELECT pq.SourceFile AS sample_id, pq.PeptideModSeq AS sequence_modified, pq.PeptideSeq AS sequence_plain, pq.PrecursorCharge AS charge, pq.RTInSecondsCenter AS rt, pq.TotalIntensity AS intensity, ps.QValue AS confidence FROM peptidequants pq NATURAL JOIN peptidescores ps'))
    n_unique_raw_files = n_distinct(p$sample_id)
    if(n_unique_raw_files == 1) {
      peptides = bind_rows(peptides, p)
      pep2prot = bind_rows(pep2prot, as_tibble(DBI::dbGetQuery(db, 'SELECT PeptideSeq AS sequence_plain, ProteinAccession AS protein_id, isDecoy AS isdecoy FROM peptidetoprotein')))
    } else {
      append_log(paste("skipping input that contains summarized data (eg. a EncyclopeDIA chromatogram library):", filename), type = "warning")
    }
    rm(p)

    # close DB connection
    DBI::dbDisconnect(db)
  }

  n_unique_raw_files = n_distinct(peptides$sample_id)
  if (length(peptides) == 0 || n_unique_raw_files < 2) {
    append_log(paste("There must be at least 2 samples in the input directory:", path), type = "error")
  }

  # convert peptide RT from seconds to minutes
  peptides$rt = peptides$rt / 60

  # discard all duplicated pep2prot entries (eg; the SQLite databases will hold distinct subsets of peptides, and thus pep2prot mappings)
  pep2prot = distinct_all(pep2prot)

  # collapse protein_id by peptide sequence. We here assume that isdecoy is a protein-level decoy annotation, thus all peptides assigned to a decoy protein are decoy entries
  proteins = pep2prot %>%
    group_by(sequence_plain) %>%
    summarise(protein_id = paste0(protein_id, collapse=";"),
              isdecoy = max(isdecoy))

  # add protein_id and decoy flag to peptides tibble
  peptides = left_join(peptides, proteins, by="sequence_plain")

  # ! re-use our generic data table parser to ensure input data standardization and maintain consistency among all data parsing routines
  # functions include; enforce data types (eg; isdecoy as boolean), format data (eg; strip extensions from filenames) and remove invalid rows (eg; intensity=0)
  # also note that, at current default settings, this function selects the 'best' peptide_id for each sequence_modified (eg; if the same sequence was found with multiple charges)
  ds = import_dataset_in_long_format(x = peptides,
                                     attributes_required = list(sample_id = "sample_id",
                                                                protein_id = "protein_id",
                                                                sequence_plain = "sequence_plain",
                                                                sequence_modified = "sequence_modified",
                                                                charge = "charge",
                                                                rt = "rt",
                                                                isdecoy = "isdecoy",
                                                                intensity = "intensity",
                                                                confidence = "confidence"),
                                     confidence_threshold = confidence_threshold,
                                     return_decoys = return_decoys,
                                     do_plot = F) # data for Cscore histograms is not available in EncyclopeDIA, so disable plotting

  ds$acquisition_mode = "dia"
  return(ds)
}
