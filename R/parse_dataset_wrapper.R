

#' A generic wrapper around data import functions included with MS-DAP
#'
#' This function simply calls the respective implementation with default parameters. In most cases it's more appropriate to explicitly call the respective import function for your data, e.g. import_dataset_skyline(). The main use-case of this function is to convenience implementation in the Graphical User Interface (GUI).
#'
#' @param path the full file path of the input file or directory
#' @param type the type of dataset you want to import
#' @export
import_dataset__generic = function(path, type) {
  if(type == "skyline_dda") {
    return(import_dataset_skyline(path, acquisition_mode = "dda"))
  }
  if(type == "skyline_dia") {
    return(import_dataset_skyline(path, acquisition_mode = "dia"))
  }
  if(type == "spectronaut") {
    return(import_dataset_spectronaut(path))
  }
  if(type == "openswath") {
    return(import_dataset_openswath(path))
  }
  if(type == "diann") {
    return(import_dataset_diann(path))
  }
  if(type == "fragpipe") {
    return(import_dataset_fragpipe(path))
  }
  if(type == "peaks") {
    return(import_dataset_peaks(path))
  }
  if(type == "openms") {
    return(import_dataset_openms_mztab(path))
  }
  if(type == "pd") {
    return(import_dataset_proteomediscoverer_txt(path))
  }
  if(type == "maxquant") {
    return(import_dataset_maxquant_evidencetxt(path))
  }
  if(type == "metamorpheus") {
    return(import_dataset_metamorpheus(path))
  }
  if(type == "encyclopedia") {
    return(import_dataset_encyclopedia(path))
  }

  append_log("unsupported dataset type", type="error")
}
