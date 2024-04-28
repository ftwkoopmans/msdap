
#' return msdap package version as a string
#'
#' simple wrapper around utils::packageVersion()
#' @export
msdap_version = function() {
  as.character(packageVersion("msdap"))
}



#' add short prettyprint labels to a protein tibble using gene symbols, if available
#'
#' @param tib tibble with at least a protein_id column, but should also contain a gene_symbols_or_id column for truly pretty-print names
#' @param maxlength max string length for each output value
#' @param shorten_ambiguous replace all ambiguous entries with the first value followed by an asterix. Example; `proteinA;proteinB` becomes  `proteinA*`
#' @return input tibble with resulting values in the label column (overwritten if this column already existed in input)
add_protein_prettyprint_label = function(tib, maxlength = 12, shorten_ambiguous = TRUE) {
  if(length(tib) == 0) {
    return(NULL)
  }
  stopifnot(is.data.frame(tib) && length(tib) > 0 && "protein_id" %in% colnames(tib))

  # default label; protein_id
  tib = tib %>% mutate(label = protein_id)

  # use gene symbols if available
  if("gene_symbols_or_id" %in% colnames(tib)) {
    tib = tib %>% mutate(label = gene_symbols_or_id)
    rows_nolabel = is.na(tib$label) | nchar(tib$label) < 2
    tib$label[rows_nolabel] = tib$protein_id[rows_nolabel]
  }

  # reduce the string length of very long labels
  rows = nchar(tib$label) > maxlength
  tib$label[rows] = paste0(substr(tib$label[rows], 1, maxlength-3), "...")

  # find labels containing semi-colons, indicating ambiguous IDs
  if(shorten_ambiguous) {
    rows = grepl(";", tib$label, fixed = T)
    if(any(rows)) {
      # strip ambiguous IDs
      tib$label = gsub(";.*", "", tib$label)
      # mark ambiguous IDs with an asterix
      tib$label[rows] = paste0(tib$label[rows], "*")
    }
  }

  return(tib)
}



#' To facilitate multiprocessing, create a Parallel Socket Cluster object
#'
#' You only need this if a) not using the quickstart function and b) running DEA algorithms msqrob or msqrobsum
#'
#' @param n_thread maximum number of threads
#' @importFrom parallel detectCores makePSOCKcluster
#' @importFrom doParallel registerDoParallel
#' @export
initialize_multiprocessing = function(n_thread = 0) {
  # if no count is specified, use number of available cores minus one
  n_max = parallel::detectCores()
  n_thread = ifelse(length(n_thread) != 1 || !is.finite(n_thread) || n_thread < 1, max(c(1, n_max-1)), ceiling(n_thread))
  # limit to core count
  n_thread = min(n_max, n_thread)
  # create a psock cluster. Parallels package warns against using the "FORK" cluster with GUI front-ends or multi-threaded libraries
  cl = parallel::makePSOCKcluster(n_thread, setup_strategy = "sequential") # the latter argument deals with a bug in parallel in R 4.0; https://github.com/rstudio/rstudio/issues/6692
  doParallel::registerDoParallel(cl)
  append_log(sprintf("using %d threads for multiprocessing", n_thread), type = "info")
  return(cl)
}



#' remove NA, non-character and characters shorter than N from array
#'
#' @param x character array
#' @param minchar integer, minimum character of elements in x
remove_by_charlength = function(x, minchar) {
  x[!is.na(x) & is.character(x) & nchar(x) >= minchar]
}



#' split array into a list of N-sized chunks
#'
#' @param x array to be chunked
#' @param chunk_size integer, maximum chunk size
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



#' CoV computation for an array of non-log values
#'
#' @param nonlog_array numeric array of non-log values
coefficient_of_variation = function(nonlog_array) {
  sqrt(expm1(sd(log(nonlog_array), na.rm = T)^2))
}



#' fast CoV computation for all rows in a matrix of natural-log transformed values
#'
#' @param matrix_with_naturallog_intensities a matrix with values transformed into natural log
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



#' convert a wide-format data.frame to a matrix
#'
#' after pivoting a long-format tibble to wide-format, use this function to easily convert to a matrix (e.g. for efficient numerical computation using the matrixStats package)
#' @examples
#' \dontrun{
#'   m = as_matrix_except_first_column(dataset$peptides %>%
#'     pivot_wider(id_cols = "peptide_id", names_from = "sample_id", values_from = "intensity")
#'   )
#' }
#'
#' @param x a wide-format data.frame or tibble where the first column are assumed to represent unique row identifiers and other columns assumed 'data' columns for the output matrix
as_matrix_except_first_column = function(x) {
  if (is_tibble(x)) {
    x = as.data.frame(x, stringsAsFactors = F)
  }
  rownames(x) = x[, 1]
  as.matrix(x[, -1, drop = F])
}



#' wrapper function around tidyverse's gather to convert a matrix to long-format tibble
#'
#' NA values are removed from output
#'
#' @param mat a matrix
#' @param value_name target column name for the values
#' @param column_name target column name for the row names
#' @param row_name target column name for the column names
#' @importFrom tidyr gather
matrix_to_long = function(mat, value_name = "value", column_name = "sample", row_name = "sequence") {
  x = as_tibble(mat)
  x$rownms = rownames(mat)
  res = x %>% gather(colname, value, -rownms, na.rm = T)
  names(res)[1] = row_name
  names(res)[2] = column_name
  names(res)[3] = value_name
  return(res)
}



#' fast replacement in a data.frame; infinite values to NA
#' https://stackoverflow.com/a/57092321
#' @param df data.frame to replace values in
replace_inf_with_NA <- function(df) {
  rapply(df, function(x) replace(x, is.infinite(x), NA), classes = "numeric", how = "replace")
}



#' fast replacement in a data.frame; zero values to NA
#' https://stackoverflow.com/a/57092321
#' @param df data.frame to replace values in
replace_zero_with_NA <- function(df) {
  rapply(df, function(x) replace(x, x==0, NA), classes = "numeric", how = "replace")
}



#' fast replacement in a vector; nonfinite to target value
#' @param x vector to replace values in
#' @param val the value to replace with
replace_nonfinite = function(x, val = NA) {
  replace(x, !is.finite(x), val)
}



#' generic input validation; all provided arguments should be numeric (not array)
#'
#' @param ... arbitrary set of parameters
#' @param minval_ optionally, a finite number indicating the minimum value for each parameter
#' @param maxval_ optionally, a finite number indicating the maximum value for each parameter
check_parameter_is_numeric = function(..., minval_ = NA, maxval_ = NA) {
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
    if (is.finite(minval_) && x < minval_) {
      append_log(sprintf("function argument '%s' must be greater than %d", n, minval_), type = "error")
    }
    if (is.finite(maxval_) && x > maxval_) {
      append_log(sprintf("function argument '%s' must be smaller than %d", n, maxval_), type = "error")
    }
  }
}



#' generic input validation; all provided arguments should be boolean/logical (not array)
#'
#' @param ... arbitrary set of parameters
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



#' generic input validation; all provided arguments should be strings (not array)
#'
#' @param ... arbitrary set of parameters
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



#' generate colours analogous to ggplot's default palette
#'
#' https://stackoverflow.com/a/8197703
#'
#' @param n number of colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
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



#' convert base-log from log2 to natural log
#'
#' log_a(x) = log_b(x) / log_b(a)
#'
#' @examples
#' \dontrun{
#'   x = 1:4
#'   cbind(input = x, as_log2 = log2(x), as_ln = log2_to_ln(log2(x)),
#'   ln_then_undo = exp(log2_to_ln(log2(x))) )
#' }
#' @param x numeric values that were log2 transformed
log2_to_ln = function(x) {
  x / log2(exp(1))
}



#' threshold finite numeric values
#'
#' @param x numeric array or matrix
#' @param threshold limit finite values below this val
threshold_numerics = function(x, threshold) {
  x[is.finite(x) & x < threshold] = threshold
  return(x)
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



#' Loess fit on finite value pairs. Ignores overfit warnings, returns input data on error
#'
#' if a robust loess fit is desired, add parameter family='symmetric' which will be passed to loess() function
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param span span argument for loess() function
#' @param remove_na_yval set all output values where is.na(y) to NA. default is TRUE
#' @param min_values minimum number of values required (if less, input parameter y is returned). default is 10
#' @param ... passed to loess() function
smooth_loess_custom = function(x, y, span=.1, remove_na_yval = TRUE, min_values = 10, ...) {
  # debug; X<<-x; Y<<-y

  # subset entries where both x and y are finite, to be used for loess fit
  x_orig = x
  y_orig = y
  flag_valid = is.finite(x) & is.finite(y)
  x = x[flag_valid]
  y = y[flag_valid]

  if(length(y) >= min_values) {
    # ignore loess warnings (eg; overfitting). Especially important for systems that halt on warning (eg; options(warn = 2) ) as this occurs on many datasets where this function is used in RT QC plots
    z = tryCatch(predict(suppressWarnings(loess(y ~ x, span = span, ...)), x_orig), error=function(e) NULL)

    # ? optionally, we could fall-back and try a wider span on error. eg; add param fallback_spans=c(0.25,0.5,1); for(span in na.omit(fallback_spans)) z=tryCatch(fit loess)... break loop on success

    # if loess fit succeeded, return predicted values
    if(length(z) == length(y_orig)) {
      if(remove_na_yval) {
        z[is.na(y_orig)] = NA
      }
      return(z)
    }
  }

  # default / fallback
  return(y_orig)

  ## v1
  # if(length(y) < 10) return(y_orig)
  # z = predict(loess(y ~ x, span = span), x_orig)
  # if(remove_na_yval) {
  #   z[is.na(y_orig)] = NA
  # }
  # return(z)
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



#' format a string of R code using the formatR and styler packages
#'
#' parameter can either be a string or an array of strings (eg; lines from R history)
#' @param rcode_lines R code as an array of characters (strings)
#' @importFrom formatR tidy_source
#' @importFrom styler cache_deactivate style_text
format_r_code = function(rcode_lines) {
  rcode_formatr = formatR::tidy_source(text=rcode_lines, width=20, output = F, comment = F, blank = F)
  styler::cache_deactivate(verbose = FALSE)
  rcode_formatr_styler = styler::style_text(rcode_formatr$text.tidy, scope = "line_breaks", strict = T)
  return(as.character(rcode_formatr_styler))
}



#' subset lines of R code (provided as strings) from the last "import_dataset" function on
#'
#' in case users run the pipeline a few times consecutively (eg; testing various parameters) without resetting R history, we only want to keep the relevant portion of the code.
#' simple/transparent way is to show everything from latest dataset import function or source() on
#'
#' @param rcode_lines R code as an array of characters (strings)
subset_relevant_code_snippet_for_report = function(rcode_lines) {
  i = grep("^[^#]+=\\s*import_dataset", rcode_lines)
  # i = grep("^\\s*source\\(|\\s*=\\s*import_dataset", rcode_lines, perl = T)
  if(length(i) > 0) {
    rcode_lines = rcode_lines[tail(i,1):length(rcode_lines)]
  }
  return(rcode_lines)
}



#' remove starting or trailing substrings that ar the same for all strings in a vector
#'
#' eg; `.txt` in all filesnames, or `<experiment data><samplename>`
#'
#' @examples
#' \dontrun{ strip_common_substring(s = unique(peptides$sample_id)) }
#' \dontrun{ strip_common_substring(s = c("sample_a.txt","sample_b1.txt","sample_xyz.txt")) }
#' \dontrun{ strip_common_substring(s = c("","sample_b1.txt","sample_xyz.txt")) }
#' \dontrun{ strip_common_substring(s = c("sample.txt","sample_b1.txt","sample_xyz.txt")) }
#' \dontrun{
#'   assertthat::are_equal(
#'     strip_common_substring(c("sample_a.txt","sample_b1.txt","sample_xyz.txt")),
#'     c("a","b1","xyz")
#'   )
#' }
#'
#' @param s array of strings
#'
#' @importFrom data.table uniqueN
#' @importFrom stringr str_split_fixed str_sub
strip_common_substring = function(s) {
  s_nchar = nchar(s)
  if(any(s_nchar <= 1) || length(s) < 2) return(s)

  ## find common prefix; convert string vector to character matrix, then check which columns are identical
  # character matrix + find similarity across vector s
  m = stringr::str_split_fixed(s, "", n = max(s_nchar))
  m_uniqueN = apply(m, 2, data.table::uniqueN)

  position_first_unique = which(m_uniqueN != 1)
  if(length(position_first_unique) > 0) {
    s = stringr::str_sub(s, position_first_unique[1])
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



#' use optim() to fit a t-distribution with fixed mu=0
#' @examples
#' \dontrun{
#' # fit normal
#' x = rnorm(1000, sd = 2)
#' fit = suppressWarnings(fit_t_dist_fixed_mu(x))
#' h = hist(x, breaks = 25, freq = F)
#' curve(dt(x / fit[1], df=fit[2]) / fit[1], col = 2, add = T)
#' }
#' \dontrun{
#' # fit t-distribution
#' x = rt(1000, df = 3)
#' fit = suppressWarnings(fit_t_dist_fixed_mu(x))
#' h = hist(x, breaks = 25, freq = F)
#' curve(dt(x / fit[1], df=fit[2]) / fit[1], col = 2, add = T)
#' }
#' @param x numeric vector of at least 10 finite values
fit_t_dist_fixed_mu = function(x) {
  mydt = function(par) {
    s = par[1]
    df = par[2]
    -sum(stats::dt((x - 0)/s, df, log = TRUE) - log(s))
  }

  x = x[is.finite(x)]
  if(length(x) < 10) return(NULL)

  fit = tryCatch(
    {
      x_sd = sd(x)
      x_mad = mad(x)
      optim(par = c(x_sd*0.5, 10), fn = mydt, method = "L-BFGS-B", lower = c(x_mad*0.5, 1), upper = c(x_sd*1.5, 20))
    },
    error = function(x) NULL
  )
  if(!is.null(fit)) return(c(sigma=fit$par[1], df=fit$par[2]))
}
