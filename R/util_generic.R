# https://stackoverflow.com/questions/36338629/escaping-special-latex-characters-in-r
escapeLatexSpecials <- function(x) {
  x <- gsub("\\", "$\\backslash$", x, fixed = T)
  x <- gsub("#", "\\\\#", x)
  x <- gsub("$", "\\\\$", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("~", "\\\\~", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("^", "\\\\^", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub(">", "$>$", x)
  x <- gsub("<", "$<$", x)
  return(x)
}



#' extract documentation for all parameters of some function
#'
#' returns a named array where values are help text and names are the function parameters
#'
#' @param func function name as a string
#' @param package package name as a string
#' @importFrom tools Rd_db
#' @export
parameter_documentation = function(func=NULL, package="msdap") {
  txt = NULL

  tryCatch({
    x = help(topic = as.character(func), package = as.character(package))
    if(length(x) > 0 && "path" %in% names(x) && file.exists(x$path) && endsWith(tolower(x$path), ".md")) {
      txt = readLines(x$path)
    }

    if(length(txt) == 0) {
      x = tools::Rd_db(as.character(package))
      # print(capture.output(print(  x  ))) # debug
      if(length(x) > 0) {
        func_name = grep(sprintf("^%s(.rd|$)", func), names(x), ignore.case = T, value = T)
        if(length(func) == 1 && length(func_name) == 1) {
          txt = x[[func_name]]
        } else {
          txt = x
        }
        txt = capture.output(print(txt))
      }
    }



    if(length(txt) > 0) {
      s = sub("^\\\\item\\{ *(.+) *} *\\{ *(.+) *}$", "\\1###\\2", txt)

      params = strsplit(s[s != txt], "###", fixed = T)
      # as matrix; matrix(unlist(params), ncol=2, byrow = T, dimnames=list(NULL, c("param","help")) )

      result = unlist(lapply(params, "[[", 2), use.names = F)
      names(result) = unlist(lapply(params, "[[", 1), use.names = F)
      return(result)
    }
  }, error=function(err) { print(err) })

  append_log("could not find package documentation file", type = "warning")
}



# assumes all are valid
contrast_to_prettyprint = function(l) {
  unlist(lapply(l, function(x) {
    paste(paste(x[[1]],collapse=","),
          paste(x[[2]],collapse=","),
          collapse=" vs ")
  }), use.names = F, recursive = F)
}



# input; a list of lists (the nested lists must be of length 2)
# @return; list of error strings
contrast_definitions_validate = function(l) {
  if(length(l) == 0) {
    return()
  }

  is_contrast_valid = unlist(lapply(l, function(x) {
    is.list(x) && length(x) == 2 && length(x[[1]]) > 0 && length(x[[2]]) > 0
  }))

  # don't even bother checking for dupes if the input list is not a valid contrast list (eg; )
  if(!all(is_contrast_valid)) {
    return("in A/B testing, there must be sample groups on both left- and right-hand side of each contrast")
  }

  ## check for duplicates, considering swaps within and between groups
  # eg; WT vs KO  ==  KO vs WT
  # eg; a,b vs c,d  ==  b,a vs c,d
  l_comparable = lapply(l, function(x) {
    paste(sort(c(paste(sort(x[[1]]),collapse=","),
                 paste(sort(x[[2]]),collapse=","))),
          collapse="###")
  })

  index_dupes = which(duplicated(unlist(l_comparable, use.names=F)))
  if(length(index_dupes) > 0) {
    return(paste("this contrast is a duplicate: ", contrast_to_prettyprint(l[index_dupes])))
  }
}


# utility function for printing R code to the logfile/report
r_variable_to_string = function(x) {
  s = capture.output(print(parse(text=deparse(x))) )
  sub("^expression\\(", "", substr(s, 1, nchar(s)-1))
}




#' To facilitate multiprocessing, create a Parallel Socket Cluster object
#'
#' You only need this if a) not using the quickstart function and b) running DEA algorithms msqrob or msqrobsum
#'
#' @param n_thread number of threads
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @export
initialize_multiprocessing = function(n_thread = 0) {
  # if no count is specified, use number of available cores minus one
  n_max = parallel::detectCores()
  n_thread = ifelse(!is.finite(n_thread) || n_thread <= 0, max(c(1, n_max-1)), n_thread)
  # limit to core count
  n_thread = min(n_max, n_thread)
  # create a psock cluster. Parallels package warns against using the "FORK" cluster with GUI front-ends or multi-threaded libraries
  cl = parallel::makePSOCKcluster(n_thread, setup_strategy = "sequential") # the latter argument deals with a bug in parallel in R 4.0; https://github.com/rstudio/rstudio/issues/6692
  doParallel::registerDoParallel(cl)
  append_log(sprintf("using %d threads for multiprocessing", n_thread), type = "info")
  return(cl)
}




#' placeholder title
#' https://stackoverflow.com/a/40622857
#' @param lst todo
#' @param envir todo
unpack_list <- function(lst, envir = parent.frame()) {
  invisible(list2env(lst, envir = envir))
}



#' sort; mixedorder the gene symbols, but put 'has no gene symbol' (protein long ID is used as placeholder) on top
#' @param s todo
#'
#' @importFrom gtools mixedorder
order_gene_symbols = function(s) {
  ss = s
  ss[grepl("^[^;]+\\|", ss)] = "" # those that have a gene symbol for the 'leading' proteingroup are not placed in front
  gtools::mixedorder(ss)
} # s=c("aa","stx12","stx1","sp|kasfsd|sdfdsf","bb;sp|kasfsd|sdfdsf"); cbind(s, s[order_gene_symbols(s)])



#' placeholder title
#' @param x todo
remove_empty = function(x) {
  x[!is.na(x) & x != ""]
}



#' split array into a list of N-sized chunks
#' @param x todo
#' @param chunk_size todo
split_array = function(x, chunk_size) {
  l = list()
  i = 1
  n = length(x)
  while (i <= n) {
    j = min(n, i + chunk_size - 1)
    l[[length(l) + 1]] = x[i:j]
    i = i + chunk_size
  }
  return(l)
}



#' placeholder title
#' @param nonlog_array todo
coefficient_of_variation = function(nonlog_array) {
  sqrt(expm1(sd(log(nonlog_array), na.rm = T)^2))
}

#' fast CoV computation
#' @param matrix_with_naturallog_intensities a matrix with values transformed into natural log
#'
#' @importFrom matrixStats rowSds
coefficient_of_variation_vectorized = function(matrix_with_naturallog_intensities) {
  sqrt(expm1( matrixStats::rowSds(matrix_with_naturallog_intensities, na.rm = T)^2 ))
}


#' placeholder title
#' @param bp todo
#' @param show_N todo
#' @param show_values todo
#' @param ... todo
boxplot_add_text = function(bp, show_N = TRUE, show_values = T, ...) {
  if (show_N) {
    text(x = seq_along(bp$n), y = bp$stats[3, ], labels = paste("n=", bp$n, sep = ""), adj = c(1, -0.5), ...)
  }
  if (show_values) {
    for (i in seq_along(bp$n)) {
      text(x = i, y = bp$stats[, i], labels = sapply(bp$stats[, i], function(x) sprintf("%.2f", x)), adj = c(-0.5, -0.5), ...)
    }
  }
}





#' placeholder title
#' @param list.of.strings todo
#' @param referenceList todo
#' @param remove.na todo
match_list_elements_to_vector = function(list.of.strings, referenceList, remove.na = T) {
  ul = unlist(list.of.strings)
  i = match(ul, referenceList)
  l = relist(i, skeleton = list.of.strings)
  if (remove.na) {
    l = lapply(l, na.omit)
  }
  return(l)
}



#' placeholder title
#' @param x todo
as_matrix_except_first_column = function(x) {
  if (is_tibble(x)) {
    x = as.data.frame(x, stringsAsFactors = F)
  }
  rownames(x) = x[, 1]
  as.matrix(x[, -1, drop = F])
}



#' placeholder title
#' @param mat todo
#' @param value_name todo
#' @param column_name todo
#' @param row_name todo
matrix_to_long = function(mat, value_name = "value", column_name = "sample", row_name = "sequence") {
  x = as_tibble(mat)
  x$rownms = rownames(mat)
  res = x %>% gather(colname, value, -rownms, na.rm = T)
  names(res)[1] = row_name
  names(res)[2] = column_name
  names(res)[3] = value_name
  return(res)
}



#' placeholder title
#' @param mat todo
#' @param groups todo
#' @param FUN todo
#' @param ... todo
aggregate_by_group = function(mat, groups, FUN, ...) {
  ugroup = unique(groups)
  res = matrix(NA, nrow(mat), length(ugroup), dimnames = list(rownames(mat), ugroup))
  for (i in seq_along(ugroup)) {
    res[, i] = apply(mat[, groups == ugroup[i], drop = F], 1, FUN, ...)
  }
  return(res)
}



#' fast replacement in a data.frame
#' https://stackoverflow.com/a/57092321
#' @param df data.frame to replace values in
replace_inf_with_NA <- function(df) {
  rapply(df, function(x) replace(x, is.infinite(x), NA), classes = "numeric", how = "replace")
}



#' fast replacement in a data.frame
#' https://stackoverflow.com/a/57092321
#' @param df data.frame to replace values in
replace_zero_with_NA <- function(df) {
  rapply(df, function(x) replace(x, x==0, NA), classes = "numeric", how = "replace")
}



#' fast replacement in a vector
#' @param x vector to replace values in
#' @param val the value to replace with
replace_nonfinite = function(x, val = NA) {
  replace(x, !is.finite(x), val)
}



#' placeholder title
#' @param ... todo
#' @param minval todo
#' @param maxval todo
check_parameter_is_numeric = function(..., minval = NA, maxval = NA) {
  varnames = lapply(substitute(list(...))[-1], deparse)
  arguments <- list(...)
  names(arguments) = unlist(varnames, use.names = F)

  for (n in names(arguments)) {
    x = arguments[[n]]
    if (length(x) != 1) {
      append_log(sprintf("function argument '%s' must be a single value", n), type = "error")
    }
    if (!is.finite(x) || !is.numeric(x)) {
      append_log(sprintf("function argument '%s' must be numeric", n), type = "error")
    }
    if (is.finite(minval) && x < minval) {
      append_log(sprintf("function argument '%s' must be greater than %d", n, minval), type = "error")
    }
    if (is.finite(maxval) && x > maxval) {
      append_log(sprintf("function argument '%s' must be smaller than %d", n, maxval), type = "error")
    }
  }
}



#' placeholder title
#' @param ... todo
check_parameter_is_boolean = function(...) {
  varnames = lapply(substitute(list(...))[-1], deparse)
  arguments <- list(...)
  names(arguments) = unlist(varnames, use.names = F)

  for (n in names(arguments)) {
    x = arguments[[n]]
    if (length(x) != 1) {
      append_log(sprintf("function argument '%s' must be a single value", n), type = "error")
    }
    if (!is.finite(x) || !is.logical(x)) {
      append_log(sprintf("function argument '%s' must be boolean", n), type = "error")
    }
  }
}



#' placeholder title
#' @param ... todo
check_parameter_is_string = function(...) {
  varnames = lapply(substitute(list(...))[-1], deparse)
  arguments <- list(...)
  names(arguments) = unlist(varnames, use.names = F)

  for (n in names(arguments)) {
    x = arguments[[n]]
    if (length(x) != 1) {
      append_log(sprintf("function argument '%s' must be a single value", n), type = "error")
    }
    if (!is.character(x)) {
      append_log(sprintf("function argument '%s' must be a string", n), type = "error")
    }
  }
}



#' placeholder title
#' @param x todo
val2ecdf = function(x) {
  q = quantile(x, probs = c(.01, .99), na.rm = T)
  x[!is.finite(x) | x < q[1]] = q[1]
  x[x > q[2]] = q[2]
  x_distr = ecdf(x)
  x_distr(x)
}



#' placeholder title
#' @param x todo
val2col = function(x) {
  y = val2ecdf(x)
  z = cut(y * 10, breaks = 0:10)
  clr = rev(RColorBrewer::brewer.pal(11, "Spectral")) # RColorBrewer::display.brewer.all()
  return(clr[as.numeric(z)])
}



#' Accept pvalues that are exactly zero as the 'best pvalue'
#'
#' @param x values to transform
#' @export
minlog10 = function(x) {
  y = -log10(x)
  y[!is.finite(x)] = NA
  y[is.finite(x) & x == 0] = max(y[is.finite(y)])
  return(y)
}



#' convert base-log from log2 to natural log. log_a(x) = log_b(x) / log_b(a)
#' example: log2(4); exp(log2_to_ln(log2(4)))
#' @param x numeric values that were log2 transformed
log2_to_ln = function(x) {
  x / log2(exp(1))
}



#' Classify an array of strings by regular expression
#'
#' For example, provide fasta headers and an array of regular expressions to find out which are yeast/human while removing ambigous entries or those matching ecoli:
#' regex = c(human="_HUMA", yeast="_YEAS", discard="_ECOL")
#'
#' @param x array of string where each regex is applied to
#' @param regex named vector of regular expressions. setting proper names is important, these are the classification labels. example: c(human="_HUMA", yeast="_YEAS", discard="_ECOL")
#' @param discard_label default label for classifications that should be discarded downstream. default: "discard"
#' @export
regex_classification = function(x, regex, discard_label = "discard") {
  stopifnot(is.vector(regex) && all(nchar(names(regex))>0))

  # allocate matrix to hold each regex' results
  mat_mask = matrix(F, nrow = length(x), ncol = length(regex), dimnames=list(NULL, names(regex)))
  # apply each regex and store in matrix
  for(i in seq_along(regex)) {
    if(!is.na(regex[i]) & regex[i] != "") {
      mat_mask[,i] = grepl(regex[i], x, ignore.case=T, perl = T) # enable perl regex so we support negative lookahead
    }
  }

  # count how many regex' match each protein / fasta header
  mat_mask_counts = rowSums(mat_mask)
  rows_single = mat_mask_counts == 1

  # default return value; discard
  classification = rep(discard_label, length(x))
  # for each regex, update classification if protein only matches this regex
  for(n in names(regex)) {
    classification[rows_single & mat_mask[,n]] = n
  }

  # proteins that match multiple regex' are discarded
  classification[mat_mask_counts > 1] = discard_label

  return(classification)
}



#' bin an array and return both the bins and their metadata (label and mean value)
#'
#' @param x numeric vector
#' @param from parameter for seq()
#' @param to parameter for seq()
#' @param length length.out parameter for seq()
make_bins = function(x, from, to, length) {
  x_mask = !is.na(x)
  b = seq(from=from, to=to, length.out = length)
  x[x_mask & x<from] = from
  x[x_mask & x>to] = to
  ct = cut(x, breaks = b, include.lowest = T) # ct[as.numeric(ct) == 1]
  # from strings that contain the range of each bin in brackets, extract both numerics and take mean value
  ct_mean_value = unlist(lapply(strsplit(gsub("[^0-9,.-]", "", levels(ct)), split = ",", fixed = T),
                                function(x) mean(as.numeric(x))))
  ct_numeric = as.numeric(ct)
  result = array(ct_numeric, dimnames = list(as.character(ct_mean_value)[ct_numeric]))

  attr(result, "bin_ids") <- seq_along(levels(ct))
  attr(result, "bin_names") <- levels(ct)
  attr(result, "bin_means") <- ct_mean_value
  return(result)
}



#' wrapper for; predict(loess(y ~ x, span = span, na.action = na.omit), x)
#'
#' also adds option to disregard predicted values if input y is NA.
#' does nothing if less than 10 values are passed.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param span argument for loess()
#' @param remove_na_yval set all output values where is.na(y) to NA
smooth_loess_custom = function(x, y, span=.1, remove_na_yval = TRUE) {
  if(length(y) < 10) return(y)
  # X<<-x; Y<<-y
  z = predict(loess(y ~ x, span = span, na.action = na.omit), x)
  if(remove_na_yval) {
    z[is.na(y)] = NA
  }
  return(z)
}



#' scale an array of numerics between 0~1, thresholding values at lower/upper quantiles
#' non-finite values are ignored (not replaced).
#' example: x=rnorm(1000); hist(x); hist(scale_between_quantiles(x, min_quantile=0.01, max_quantile = 0.99))
#'
#' edge-case; if either lower or upper quantile could not be determined, all finite values of x will be set to 1
#' edge-case; if all non-finite input values are the same, the upper and lower quantiles are the same. All finite values of x will be set to 1
#'
#' @param x values to transform
#' @param min_quantile lower quantile, used as parameter in quantile()
#' @param max_quantile upper quantile, used as parameter in quantile()
scale_between_quantiles = function(x, min_quantile = 0.01, max_quantile = 0.99) {
  # basic input validation
  stopifnot(is.numeric(x) && min_quantile>0 && max_quantile<1 && max_quantile > min_quantile)

  x_mask = is.finite(x)
  q = quantile(x[x_mask], probs = c(min_quantile, max_quantile))

  # edge-case; if all non-finite input values are the same, the upper and lower quantiles are the same. All finite values of x will be set to 1
  if(any(!is.finite(q)) || q[2] <= q[1]) {
    x[x_mask] = 1
    return(x)
  }

  x[x_mask & x < q[1]] = q[1]
  x[x_mask & x > q[2]] = q[2]
  # subtract lower threshold, then divide by upper threshold (which we adjust for the lower threshold as well)
  (x-q[1]) / (q[2] - q[1])
}



#' placeholder title
#' parameter can either be a string or an array of strings (eg; lines from R history)
#' @param rcode todo
#'
#' @importFrom formatR tidy_source
#' @importFrom styler cache_deactivate style_text
format_r_code = function(rcode) {
  rcode_formatr = formatR::tidy_source(text=rcode, width=20, output = F, comment = F, blank = F)
  styler::cache_deactivate(verbose = FALSE)
  rcode_formatr_styler = styler::style_text(rcode_formatr$text.tidy, scope = "line_breaks", strict = T)
  return(as.character(rcode_formatr_styler))
}



#' in case users run the pipeline a few times consecutively (eg; testing various parameters) without resetting R history, we only want to keep the relevant portion of the code.
#' simple/transparent way is to show everything from latest dataset import function or source() on
#'
#' @param code_lines todo
subset_relevant_code_snippet_for_report = function(code_lines) {
  i = grep("^[^#]+=\\s*import_dataset", code_lines)
  # i = grep("^\\s*source\\(|\\s*=\\s*import_dataset", code_lines, perl = T)
  if(length(i) > 0) {
    code_lines = code_lines[tail(i,1):length(code_lines)]
  }
  return(code_lines)
}



# # naive approach to printing R code cleanly for a report
# # removing comments and keeping only minimal whitespace
# # assumes 'valid' R code, zero quality control or input validation
# r_code_format_simple = function(rcode_lines, force_newlines = FALSE) {
#   if(length(rcode_lines) <= 2) return(rcode_lines)
#   # remove leading and trailing whitespace, and code comments
#   s = r_code_remove_comment_simple(gsub("(^\\s+)|(\\s+$)", "", rcode_lines))
#
#   # keep track of the identation level; difference in open- and close-parentheses indicates if some code line is supposed to be idented (eg; params within function call)
#   open_parentheses = stringr::str_count(s, stringr::fixed("("))
#   close_parentheses = stringr::str_count(s, stringr::fixed(")"))
#   ident_level_shift = open_parentheses - close_parentheses
#   ident_level = cumsum(ident_level_shift)
#   fix = c(0, head(ident_level, -1))
#   # if there is 1+ close parentheses, and nothing else, on this row; reduce identation by one level
#   s_only_close_parentheses = close_parentheses> 0 & !grepl("[^)]", s)
#   fix[fix > 0 & s_only_close_parentheses] = fix[fix > 0 & s_only_close_parentheses] - 1
#   fix[fix < 0] = 0 # should never occur, but add for robustness anyway (eg; if invalid code is provided or whatever. this implementation is simple and naive)
#   # debug; cbind(s, ident_level_shift, ident_level, fix)
#
#   if(force_newlines) {
#     # remove all newlines, then add one before each new statement
#     flag = s != ""
#     # apply identation
#     s = paste0(strrep("  ", fix[flag]), s[flag])
#     s2 = NULL
#     for(line in s) {
#       if(grepl("^[a-zA-Z]", line))
#         line = c("", line)
#       s2 = c(s2, line)
#     }
#     return(s2)
#   } else {
#     # invalidate double-newlines and newlines within function calls
#     valid = rep(T, length(s))
#     for(i in 2:(length(s)-1)) {
#       # valid; if not empty OR previous line is not empty and not idented (eg; this line is argument in function call)
#       valid[i] = s[i] != "" || (s[i-1] != "" && fix[i] == 0)
#     }
#     # apply identation
#     s = paste0(strrep("  ", fix), s)
#     return(s[valid])
#   }
# }
#
#
#
# # given some lines of R code, find comments and remove
# # takes strings into account to handle hashtag symbols used in actual code (eg; nested in some regex)
# # test line; s = c(" # test", " a = 1 # comment", 'sub("#.*", "", "string ## comment")         # comment'); cbind(s, r_code_remove_comment_simple(s))
# r_code_remove_comment_simple = function(x) {
#   x_has_hashtag = grepl("#", x, fixed=T)
#   if(!any(x_has_hashtag)) return(x) # nothing to do
#
#   s = x[x_has_hashtag]
#   # regex to remove quoted parts by strings of equal length (handling both single and double quotes)
#   s = stringr::str_replace_all(s, '"[^"]*"', function(string_to_replace) strrep(" ", nchar(string_to_replace)))
#   s = stringr::str_replace_all(s, "'[^']*'", function(string_to_replace) strrep(" ", nchar(string_to_replace)))
#   # identify first comment char in updated string
#   i = stringr::str_locate(s, "\\s*#")[,1]
#   if(!any(!is.na(i))) return(x) # nothing to do
#   # those where index i is not NA, strip from input string
#   x[x_has_hashtag][!is.na(i)] = stringr::str_sub(x[x_has_hashtag][!is.na(i)], 1, na.omit(i)-1)
#   return(x)
# }



#' remove starting or trailing substrings that ar the same for all strings in a vector. eg; .txt in all filesnames, or <experiment data><samplename>
#' strip_common_substring(s = unique(peptides$sample_id))
#' strip_common_substring(s = c("sample_a.txt","sample_b1.txt","sample_xyz.txt"))
#' strip_common_substring(s = c("","sample_b1.txt","sample_xyz.txt"))
#' strip_common_substring(s = c("sample.txt","sample_b1.txt","sample_xyz.txt"))
#' assertthat::are_equal(strip_common_substring(c("sample_a.txt","sample_b1.txt","sample_xyz.txt"), c("a","b1","xyz")))
#'
#' @param s todo
#'
#' @importFrom data.table uniqueN
#' @importFrom stringr str_split_fixed str_sub
strip_common_substring = function(s) {
  s_nchar = nchar(s)
  if(any(s_nchar <= 1)) return(s)

  ## find common prefix; convert string vector to character matrix, then check which columns are identical
  # character matrix + find similarity across vector s
  m = stringr::str_split_fixed(s, "", n = max(s_nchar))
  m_uniqueN = apply(m, 2, data.table::uniqueN)

  position_first_unique = which(m_uniqueN != 1)[1]
  if(length(position_first_unique) == 1) {
    s = stringr::str_sub(s, position_first_unique)
  }


  ## analogously for trailing, take last N characters from s (after prefix filtering) and use the same 'identity in character matrix' trick
  # take trailing
  s_nchar = nchar(s)
  n_trail = min(s_nchar)
  s_end = stringr::str_sub(s, -n_trail) # s_end = stringr::str_sub(s, start = s_nchar - n_trail + 1, end=s_nchar)
  # character matrix + find similarity across vector s
  m_end = stringr::str_split_fixed(s_end, "", n = n_trail)
  m_end_uniqueN = apply(m_end, 2, data.table::uniqueN)
  mask_end = m_end_uniqueN == 1
  if(any(mask_end)) {
    position_last_unique = max(0, tail(which(mask_end != 1), 1))
    n_remove_end = n_trail - position_last_unique
    s = stringr::str_sub(s, 1, -n_remove_end - 1)
  }

  return(s)
}
