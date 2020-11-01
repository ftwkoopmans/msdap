
#' placeholder title
#' @param f todo
remove_file_if_exists = function(f) {
  if(file.exists(f)) {
    if(!suppressWarnings(file.remove(f))) {
      append_log(sprintf("File already exists, and unable to overwrite: '%s' -->> is it currently opened?", f), type = "error")
    }
  }
}



#' placeholder title
#' @param path todo
#' @param file todo
path_exists = function(path, file) {
  f = paste0(path, "/", file)
  if (!file.exists(f)) {
    append_log(sprintf('Cannot find file "%s" in directory "%s"', file, path), type = "error")
  }
  return(f)
}



#' purposely naive approach to keep things simple for our internal usecases; strip a . followed by 0-5 letters from end of the string
#' @param f f can either be a filename or a full path, and optionally carry a file extension. eg; remove_file_extension_from_path("file.txt"); remove_file_extension_from_path("C:/temp/file")
remove_file_extension_from_path = function(f) {
  sub("\\.[a-zA-Z]{0,5}$", "", f)
}



#' placeholder title
#' @param headers todo
#' @param l list of mappings, where list names are target names and values are arrays of patterns that should match provided table headers
#' @param error_on_missing todo
#' @param allow_multiple todo
map_headers = function(headers, l, error_on_missing = FALSE, allow_multiple = FALSE) {
  # convert to lowercase and apply regex to remove non-alphanumeric characters; makes matching more robust in most use-cases
  headers_clean = gsub("[^[:alnum:]]", "", tolower(headers))
  l_clean = lapply(l, function(x) gsub("[^[:alnum:]]", "", tolower(x)))
  # map attribute list to input and remove non-matches
  l_match = lapply(l_clean, match, headers_clean)
  l_match = lapply(l_match, na.omit)

  if(error_on_missing && any(lengths(l_match) == 0)) {
    append_log(paste('No matches for required columns;', paste(names(l)[lengths(l_match) == 0], collapse = ", ")), type = "error")
  }
  if(!allow_multiple && any(lengths(l_match) > 1)) {
    append_log(paste('Multiple items match for columns;', paste(names(l)[lengths(l_match) > 1], collapse = ", ")), type = "error")
  }

  # automatically removes elements that have no matches
  unlist(lapply(l_match, head, 1), use.names = T)
  # l_match_as_array = unlist(lapply(l_match, function(x) ifelse(length(x)==0, NA, x[1])))
  # return(list(indices = l_match_as_array, flag_missing = lengths(l_match) == 0, flag_multiple = lengths(l_match) > 1))
}


#' combine PDF files
#'
#' the function pdf_combine() from pdftools/qpdf package has issues (on windows) when combining hundreds of plots. error: "Too many open files"
#' so we here create a wrapper where we combine source files in batches
#' some references; https://stackoverflow.com/questions/40810704/error-with-r-function-download-file-too-many-open-files
#' @param input todo
#' @param output todo
#'
#' @importFrom pdftools pdf_combine
pdf_combine_chunks = function(input, output) {
  # chunks of at most 200 files  &  respective temporary filenames (re-using first input file to get unique file / for unique purpose)
  input_chunks = split(input, ceiling(seq_along(input) / 200))
  output_temp = paste0(input[1], ".chunk", seq_along(input_chunks))
  # combine files within each chunk
  for(i in seq_along(input_chunks)) {
    pdftools::pdf_combine(input_chunks[[i]], output_temp[i])
  }
  # finally, combine chunks and remove temp files
  pdftools::pdf_combine(output_temp, output)
  res = file.remove(output_temp)
}
