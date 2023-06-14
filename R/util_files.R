
#' append file to path, then call file_check() for validation of the file. Does not validate illegal characters in the path
#'
#' returns a (potentially cleaned) file path
#'
#' @param path a directory on this computer. eg; "C:/temp" on windows or "/home/user1" on unix
#' @param file filename that should be appended to the path
#' @param strict boolean indicating 'strict mode' for checking the filename, if TRUE; only allow alphanumeric, underscore and '-'
path_append_and_check = function(path, file, strict = FALSE) {
  file_check(paste(c(path, file), collapse = "/"), strict = strict)
}



#' strip illegal characters from filename, then check for maximum path length
#'
#' If longer than N characters (Windows filesystem max path length), try with md5 hashed filename
#'
#' note;
#' in contrast to using a timestamp or substring+index, this solution yields stable filenames. Should never fail if provided path/dir is reasonable
#'
#' @examples
#' \dontrun{
#'   x = c("C:/temp/file*.txt", "test=2", "test=2.abc",
#'          paste0("c:/temp/",paste(rep("a",300),collapse=""),".txt"))
#'   cbind(x, file_check(x))
#' }
#' \dontrun{
#'   file_check("C:/temp/") # error, trailing slashes not allowed
#' }
#' @param file filename to check (may include a path)
#' @param strict boolean indicating 'strict mode' for checking the filename, if TRUE; only allow alphanumeric, underscore and '-'
#' @importFrom openssl md5
file_check = function(file, strict = FALSE) {
  file_asis = file
  file = path_clean_slashes(file)

  # not a file if it ends with a /
  file_error = grepl("/$", file)
  if(any(file_error)) {
    append_log(paste0("not valid filenames (trailing / not allowed);\n", paste(file[file_error], collapse = "\n")), type="error")
  }

  file_fail = nchar(file) > 260 # Windows path length limit
  file[!file_fail] = filename_strip_illegal_characters(file[!file_fail], strict = strict)

  # attempt to create a shorter path using md5 hashed filename (32 chars @ md5, perhaps shorter than intended filename. eg; auto generated filenames for statistical contrasts)
  for(i in which(file_fail)) {
    # note; cannot use dirname on rediculously long paths
    # dir is everything up to, and including, first forward slash  (note that earlier we already replaced all slashes with /)
    f_dir = gsub("/[^/]+$", "/", file[i])
    f_dir[f_dir == file[i]] = "" # no regex match = no slash in current file path = empty dirname

    # only if the directory length + expected length of md5 hash and file extension doesn't exceed limit
    if(nchar(f_dir) + 32 + 5 <= 260) {
      # directory including trailing slash (if any)  +  md5 hash  +  extension is still missing)
      file[i] = paste0(f_dir, openssl::md5(file[i]), filename_extract_extension(file[i]))
      append_log(sprintf('output file has a total path length longer than 260 characters; "%s". This was shortened to "%s"', file_asis[i], file[i]), type="warning")
    } else {
      append_log(paste("failed to create output file, total path length longer than 260 characters. Choose an output_dir with shorter path length, eg; C:/data/myproject (Windows) or /home/user/myproject (unix). Attempted output file (longer than 260 characters, so invalid);", file_asis[i]), type="error")
    }
  }

  return(file)
}



#' remove file if exists
#'
#' @param f full path to some file
remove_file_if_exists = function(f) {
  if(file.exists(f)) {
    if(!suppressWarnings(file.remove(f))) {
      append_log(sprintf("File already exists, and unable to overwrite: '%s' -->> is it currently opened?", f), type = "error")
    }
  }
}



#' Combine path and filename, then check if filename exists. Optionally, checks for present of compressed variant of this filename
#'
#' @param path a directory on this computer. eg; "C:/temp" on windows or "/home/user1" on unix. Can also provide full path here and leave 'file' param as NULL
#' @param file filename within the 'path'
#' @param try_compressed if the requested filename is not available, try again with added file extensions known for compression (default: FALSE). See `read_textfile_compressed()` function docs for supported extensions
#' @param silent on error, return an empty string instead of throwing an error (default: FALSE)
path_exists = function(path, file = NULL, try_compressed = FALSE, silent = FALSE) {
  stopifnot(length(path) == 1 && length(file) <= 1)
  f = paste(c(path, file), collapse = "/") # f = full path to test

  if(file.exists(f)) { # if found, return f
    return(f)
  }

  if(try_compressed) {
    query_filename = basename(f)
    # directory to search files in
    query_dir = dirname(f)
    if(query_dir == f) {
      query_dir = getwd()
    }

    # candidate files in target directory
    candidates = data.frame(files = dir(path = query_dir, full.names = T, recursive = F, ignore.case = T, include.dirs = F), stringsAsFactors = F) %>%
      mutate(
        filenames = basename(files),
        # regex search for filenames with supported compression extensions appended
        filenames_noext_archive = sub("\\.(zip|gz|bz2|xz|7z|zst|lz4)$", "", filenames, ignore.case = T)
      ) %>%
      # only files that actually have compression file extension; if regex remove of extension failed, not a candidate file
      filter(filenames != filenames_noext_archive)

    if(nrow(candidates) > 0) {
      # match user's query file
      if(query_filename %in% candidates$filenames_noext_archive) {
        return(candidates$files[match(query_filename, candidates$filenames_noext_archive)[1]])
      }

      # analogous, but try with the original file extension removed as some compression tools situationally shorten output filename
      # e.g. C:/temp/test.tsv -->> ZIP -->> C:/temp/test.zip  (instead of C:/temp/test.tsv.zip)
      query_filename_noext = remove_file_extension_from_path(query_filename)
      if(query_filename_noext %in% candidates$filenames_noext_archive) {
        return(candidates$files[match(query_filename_noext, candidates$filenames_noext_archive)[1]])
      }
    }
  }

  # fallthrough; failed to find file
  if(silent) {
    return("")
  } else {
    append_log(sprintf('Cannot find file "%s"', f), type = "error")
  }
}



#' remove a single file extension from a path or filename
#'
#' purposely naive approach to keep things simple for our internal usecases; strip a . followed by 0-5 characters from end of the string
#'
#' @param f f can either be a filename or a full path, and optionally carry a file extension. eg; remove_file_extension_from_path("file.txt"); remove_file_extension_from_path("C:/temp/file")
remove_file_extension_from_path = function(f) {
  sub("\\.[a-zA-Z0-9]{0,5}$", "", f) # use sub, not gsub, to match exactly once
}



#' returns a regex string that can be used to remove extensions from mass-spec raw files
#'
regex_rawfile_strip_extension = function() {
  # update; greedy replace of extensions (e.g. ".wiff.dia" -> "")
  return("(.*(\\\\|/))|((\\.(mzML|mzXML|WIFF|RAW|htrms|dia|d|zip|gz|bz2|xz|7z|zst|lz4))+$)")
  # return("(.*(\\\\|/))|(\\.(mzML|mzXML|WIFF|RAW|htrms|dia|d|zip|gz|bz2|xz|7z|zst|lz4)$)")
}



#' extract the trailing file extension from a file or path
#'
#' purposely naive approach to keep things simple for our internal usecases; strip a . followed by 1-5 characters from end of the string
#'
#' @param f f can either be a filename or a full path, and optionally carry a file extension. eg; remove_file_extension_from_path("file.txt"); remove_file_extension_from_path("C:/temp/file")
#' @return extension including the '.' or empty string if not found
filename_extract_extension = function(f) {
  ext = sub(".*(\\.[a-zA-Z0-9]{1,5})$", "\\1", f) # use sub, not gsub, to match exactly once
  ifelse(ext != f, ext, "")
}



#' regex replace back slashes and redundant forward/back-slashes with a single forward slash
#'
#' @param f full path to some file
path_clean_slashes = function(f) {
  gsub("[/\\]+", "/", f)
}



#' from filename minus extension, replace illegal characters (allowed are alphanumeric, underscore and dash)
#'
#' @examples
#' \dontrun{
#'   # nothing to do
#'   filename_strip_illegal_characters("test.txt")
#'   filename_strip_illegal_characters("C:/temp/test.txt")
#'   # fix filename
#'   filename_strip_illegal_characters("C:/temp/test !~: file.txt", strict = F)
#'   # example emphasizes that this function only applies to filename
#'   filename_strip_illegal_characters("C:/temp/this is ok/not ok.txt", strict = T)
#' }
#' @param f f can either be a filename or a full path
#' @param strict boolean, if TRUE; only allow alphanumeric, underscore and '-'
#' @param replacement string replacement value for removed characters
filename_strip_illegal_characters = function(f, strict = FALSE, replacement = "_") {
  f_file = basename(f) # actual filename, in case f is a path
  f_dir = dirname(f)
  f_dir_none = f_dir == "." | f_dir == ""
  f_dir[f_dir_none] = ""
  f_dir[!f_dir_none] = paste0(f_dir[!f_dir_none], "/")

  f_file_new = f_file
  if(strict) {
    f_file_new = paste0(gsub("[^a-zA-Z0-9_-]", replacement, remove_file_extension_from_path(f_file)), filename_extract_extension(f_file))
  } else {
    f_file_new = paste0(gsub("[^a-zA-Z0-9 .;#_-]", replacement, remove_file_extension_from_path(f_file)), filename_extract_extension(f_file))
  }

  paste0(f_dir, f_file_new)
}



#' Read text files as vector of lines or a table, supports compressed files
#'
#' Supported compression formats; .zip|.gz|.bz2|.xz|.7z|.zst|.lz4
#' Using this function with other compression formats, or any of these formats but with unlisted file extensions, will not work
#'
#' note; the archive R package has some bugs still so we try to use the readr package for common compression formats
#' e.g. bug while reading ZIP files; "Error: The size of the connection buffer (131072) was not large enough". Related bugreport @ https://github.com/tidyverse/vroom/issues/361
#'
#' note; another bug is with reading Zstd archives compressed with --long
#' decompressing with current implementation results in;
#' Error: archive_read.cpp:102 archive_read_open1(): Zstd decompression failed: Frame requires too much memory for decoding
#' Related bugreport @ https://github.com/libarchive/libarchive/issues/1795
#'
#' @param file path to input file
#' @param as_table boolean, TRUE = results from data.table::fread(), FALSE = results from readr::read_lines() (default)
#' @param skip_empty_rows ignore empty rows (default: FALSE)
#' @param nrow optionally, integer specifying to only read first N rows
#' @param ... sent to data.table::fread()
#' @return NULL if path doesn't exist or file could not be read / decompressed (warnings/errors are silent)
#' @importFrom archive file_read archive_read
#' @importFrom readr read_lines
#' @importFrom data.table fread
#' @export
read_textfile_compressed = function(file, as_table = FALSE, skip_empty_rows = FALSE, nrow = -1, ...) {
  # TODO: generate warning when as_table && skip_empty_rows -->> disabled
  if(length(file) != 1 || !is.character(file) || !file.exists(file)) {
    return(NULL)
  }
  success = FALSE
  x = con = NULL
  # file extensions
  file_len = nchar(file)
  ext2 = tolower(substring(file, file_len-2, file_len))
  ext3 = tolower(substring(file, file_len-3, file_len))
  ext_use_archive = any(c(ext2, ext3) %in% c(".7z",".lz4", ".zst")) # extensions we need to use libarchive
  ext_use_readr = any(c(ext2, ext3) %in% c(".gz",".xz", ".zip", ".bz2")) # supported by readr::read_lines()
  ext_assume_compressed = ext_use_archive || ext_use_readr


  suppressWarnings(tryCatch({
    ### read lines
    # if not a table, or we need a table but the file is compressed; x = read lines @ file
    if(!as_table || ext_assume_compressed) {
      if(!ext_use_archive) {
        x = readr::read_lines(file, n_max = nrow, skip_empty_rows = skip_empty_rows, progress = FALSE)
      } else {
        ### use the archive package to create a connection to the compressed file
        # for some reason, archive::file_read doesn't work for 7z files but archive::archive_read does
        if(ext2 == ".7z") {
          con = archive::archive_read(file)
        } else {
          con = archive::file_read(file, mode = "r", filter = NULL, options = character())
        }
        # read libarchive text connection. seems less error prone when setting lazy=F and no threading when using libarchive connections
        tryCatch({
          x = readr::read_lines(con, n_max = nrow, skip_empty_rows = skip_empty_rows, progress = FALSE, num_threads = 1, lazy = F)
        }, error = function(e) {
          append_log(paste0("failed to read file in read_textfile_compressed():\n", conditionMessage(e)), type = "error")
        })
      }
    }

    ### user requested a table
    if(as_table) {
      if(!ext_assume_compressed) {
        x = data.table::fread(file = file, nrows = nrow, ...) #, blank.lines.skip = FALSE
      }
      if(ext_assume_compressed && length(x) > 0) {
        # bugfix; data.table::fread() doesn't accept input that is a single line without end-of-line character
        if(length(x) == 1) {
          x = paste0(x, "\n")
        }
        # input file is compressed, assume we already uncompressed and read all lines
        x = data.table::fread(text = x, ...) # don't need blank skip nor nrow argument, already done @ read_lines()
      }
    }

    success = TRUE
  },
  error = function(e) {
    append_log(paste0("failed to read file in read_textfile_compressed():\n", conditionMessage(e)), type = "error")
  },
  # read_line should've closed the file connection, but check anyway  (trycatch to deal with closed or NULL connections)
  finally = tryCatch({if(length(con)!=0) close(con)}, error = function(e) NULL)
  ))

  if(success == TRUE) {
    return(x)
  }
}



#' robust matching of a list of desired columns to a character array representing column headers
#'
#' instead of literal string matching between 'headers' and values in l,
#' strip/cleanup both by removing all non-alphanumeric characters and convert to lowercase.
#'
#' @param headers character array representing column names
#' @param l list of mappings, where list names are target column names and values are character arrays that are (ordered) valid values to represent a column
#' @param error_on_missing throw error if any element in l fails to match
#' @param allow_multiple throw error if any element in l has multiple matches
map_headers = function(headers, l, error_on_missing = FALSE, allow_multiple = FALSE) {
  if(length(l) == 0) {
    return()
  }
  # convert to lowercase and apply regex to remove non-alphanumeric characters; makes matching more robust in most use-cases
  headers_clean = gsub("[^[:alnum:]]", "", tolower(headers))
  l_clean = lapply(l, function(x) gsub("[^[:alnum:]]", "", tolower(x)))
  # map attribute list to input and remove non-matches
  l_match = lapply(l_clean, match, headers_clean)
  l_match = lapply(l_match, na.omit)

  log_quote_fields = function(x) paste(sapply(x, function(y) paste0('"',y,'"'), simplify = T), collapse = ", ")
  if(error_on_missing && any(lengths(l_match) == 0)) {
    err = NULL
    for(i in which(lengths(l_match) == 0)) {
      err = c(err, sprintf('field: "%s" -> could not find respective column names: %s', names(l)[i], log_quote_fields(l[[i]]) ))
    }
    append_log(paste0('Error while parsing headers of input table;\n', err), type = "error")
    # append_log(paste('No matches for required columns;', paste(names(l)[lengths(l_match) == 0], collapse = ", ")), type = "error")
  }
  if(!allow_multiple && any(lengths(l_match) > 1)) {
    err = NULL
    for(i in which(lengths(l_match) > 1)) {
      err = c(err, sprintf('field: "%s" -> searched for column names: %s -> found: %s', names(l)[i], log_quote_fields(l[[i]]), log_quote_fields(l_match[[i]]) ))
    }
    append_log(paste0('Ambiguous matching of column names in the input table;\n', err), type = "error")
    # append_log(paste('Multiple items match for columns;', paste(names(l)[lengths(l_match) > 1], collapse = ", ")), type = "error")
  }

  # automatically removes elements that have no matches
  unlist(lapply(l_match, head, 1), use.names = T)
  # l_match_as_array = unlist(lapply(l_match, function(x) ifelse(length(x)==0, NA, x[1])))
  # return(list(indices = l_match_as_array, flag_missing = lengths(l_match) == 0, flag_multiple = lengths(l_match) > 1))
}



#' Efficiently read a table based on some specification of expected headers
#'
#' Uses map_headers() to parse and map lists of expected column names and their designated names
#'
#' @param file full file path for a CSV/TSV data table
#' @param attributes_required see `map_headers()`
#' @param attributes_optional see `map_headers()`
#' @param regex_headers regex to be applied to all columns; these will be included
#' @param as_tibble_type return data.table::fread() results as a tibble (DEFAULT), instead of data.table format
read_table_by_header_spec = function(file, attributes_required, attributes_optional = NULL, regex_headers = NULL, as_tibble_type = TRUE) {
  # read first line from file; read 1 line into data.frame without colnames with first row having all values, then convert to character array
  # this leverages data.table to infer what character was used to delimit columns
  headers = as.character( read_textfile_compressed(file, as_table = T, nrow = 1, header = F, data.table = F) )
  ## analogous, but requires specification of the delimiter used;
  # headers = unlist(strsplit(read_textfile_compressed(file, as_table = F, nrow = 1), delim), recursive = F, use.names = T)

  # apply header mapping function
  map_required = map_headers(headers, attributes_required, error_on_missing = T, allow_multiple = T)
  map_optional = map_headers(headers, attributes_optional, error_on_missing = F, allow_multiple = T)

  # if a regex for columns is provided, find them in the headers array
  map_regex = NULL
  if(length(regex_headers) > 0) {
    map_regex = grep(regex_headers, headers, ignore.case = T)
    names(map_regex) = headers[map_regex]
    if(length(map_regex) == 0) {
      append_log(sprintf('Cannot find columns matching "%s" in file "%s" input file', regex_headers, file), type = "error")
    }
  }

  # collect all column indices and their desired names
  col_indices = c(map_required, map_optional, map_regex)

  # we don't allow duplicate matches
  col_indice_dupe_names = names(col_indices)[col_indices %in% col_indices[duplicated(col_indices)]]
  if(length(col_indice_dupe_names) > 0) {
    append_log(sprintf('"%s" are duplicate column names in the parsing specification applied to file "%s"', paste(col_indice_dupe_names, collapse = ", "), file), type = "error")
  }

  # only read columns of interest to speed up file parsing
  # don't parse first row as header and then skip it to avoid issues with tables that have mismatched headers (e.g. MetaMorpheus files that have extra trailing tab on header row)
  result = read_textfile_compressed(file, skip_empty_rows = F, as_table = T, select = as.integer(col_indices), header = F, skip = 1, stringsAsFactors = F, data.table = !as_tibble_type)
  colnames(result) = names(col_indices) # overwrite column names from file with the desired names from column specification
  if(as_tibble_type) {
    result = as_tibble(result)
  }

  return(result)
}



#' combine PDF files
#'
#' the function pdf_combine() from pdftools/qpdf package has issues (on windows) when combining hundreds of plots. error: "Too many open files"
#' so we here create a wrapper where we combine source files in batches
#' some references; https://stackoverflow.com/questions/40810704/error-with-r-function-download-file-too-many-open-files
#'
#' @param input array of paths to PDF files
#' @param output path to desired output PDF file
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
