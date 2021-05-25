
#' Guess experiment data from a set of filenames
#'
#' Assuming input files are somewhat well formatted,
#' we here try to detect sample groups and other metadata to populate metadata table for user's convenience
#'
#'
#' example; wt1, wt2, ko1, ko2, ...
#' example; wt_1, wt_2, ko_1, ko_2, ...
#' example; control_day1_1, control_day1_2, cond1_day1_1, cond1_day2_2, cond1_day2_3, ...
#'
#' @importFrom stringr str_count str_split_fixed
metadata_matrix_from_filenames = function(x) {
  x = tolower(x)

  # define some common separation characters, then count how often each occurs. Most common is the main sep char, others are alternatives downstream
  split_chars = count_sep_char(x)

  if(nrow(split_chars) == 0) {
    m = as.matrix(x)
  } else {
    m = stringr::str_split_fixed(x, split_chars$char_regex[1], n=Inf)
    m = shift_empty(m)
  }

  m = split_recursively(m, alt_sep_chars = split_chars$char_regex[-1])
  m = shift_empty(m)
  m = m[,apply(m, 2, n_distinct) > 1,drop=F] # remove columns without unique values
  # recognize dates formatted as 8 consecutive integers and rewrite them as "-" delimited
  m = apply(m, 2, reformat_date)

  colnames(m) = paste0("x",1:ncol(m))
  return(m)

  ## various test data;
  # x = c("KO1","KO2","KO3","WT01","WT02","WT03")
  # x = c("KO01","KO02","KO03","KO04","KO05","KO06","WT01","WT02","WT03","WT04","WT4test","WT05","WT06")
  # x = c("1s_2min","1s_2plus","1s_5min","1s_5plus","1s_8plus","1s_9min","1s_11min","1s_11plus","1s_14min","1s_14plus","3s_3min","3s_3plus","3s_6min","3s_6plus","3s_9min","3s_9plus","3s_12min","3s_12plus","3s_15min","3s_15plus","ns_4min","ns_4plus","ns_7min","ns_7plus","ns_10min","ns_10plus","ns_13min","ns_13plus")
  # x = c("WT_HC1","AD_NS1","AD_DS1","AD_DS2","WT_HC2","AD_HC1","AD_HC2","WT_NS2","AD_NS2","AD_NS3","WT_DS1","WT_DS2","AD_DS3","WT_HC3","WT_HC4","AD_HC3","AD_HC4","WT_NS3","WT_NS4","AD_NS4","AD_NS5","WT_DS3","WT_DS4","AD_DS4","WT_HC5","WT_HC6","AD_HC5","AD_HC6","WT_NS5","WT_NS6","AD_NS6","AD_NS7","WT_HC7","WT_HC8","AD_HC7","AD_HC8","WT_NS7","WT_NS8","AD_DS5","AD_DS6","WT_DS5","AD_NS8","WT_DS6","WT_DS7","WT_NS1")
  # x = c("18_15_bioni_neuron_06","18_15_bioni_neuron_glia_03","18_15_C001_neuron_07","18_15_C001_neuron_glia_04","18_15_glia_glia_01","18_15_GM_neuron_05","18_15_GM_neuron_glia_02","18_42_bioni_neuron_13","18_42_bioni_neuron_glia_10","18_42_C001_neuron_14","18_42_C001_neuron_glia_11","18_42_C002_neuron_15","18_42_glia_glia_08","18_42_GM_neuron_12","18_42_GM_neuron_glia_09","19_15_bioni_neuron_21","19_15_bioni_neuron_glia_18","19_15_C001_neuron_22","19_15_C001_neuron_glia_19","19_15_glia_glia_16","19_15_GM_neuron_20","19_15_GM_neuron_glia_17","19_42_bioni_neuron_28","19_42_bioni_neuron_glia_25","19_42_C001_neuron_glia_26","19_42_C002_neuron_29","19_42_glia_glia_23","19_42_GM_neuron_27","19_42_GM_neuron_glia_24","20_15_bioni_neuron_35","20_15_bioni_neuron_glia_32","20_15_C001_neuron_36","20_15_C001_neuron_glia_33","20_15_glia_glia_30","20_15_GM_neuron_34","20_15_GM_neuron_glia_31","20_42_bioni_neuron_42","20_42_bioni_neuron_glia_39","20_42_C001_neuron_43","20_42_C001_neuron_glia_40","20_42_C002_neuron_44","20_42_glia_glia_37","20_42_GM_neuron_41","20_42_GM_neuron_glia_38")
  # x = c("glia_C1_arrow_01","C1_arrow_25","glia_C1_cargo_02","C1_cargo_26","glia_C1_dream_03","C1_dream_27","glia_C2_dance_04","C2_dance_28","glia_C2_fire_05","C2_fire_29","glia_C2_gold_06","C2_gold_30","glia_C3_alpha_07","C3_alpha_31","glia_C3_etch_08","C3_etch_32","glia_C3_hint_09","C3_hint_33","glia_P3_carrot_10","P3_carrot_34","glia_P3_force_11","P3_force_35","P3_holly_12","glia_P3_holly_36","glia_P4_axis_13","P4_axis_37","glia_P4_cold_14","P4_cold_38","glia_P4_dawn_15","P4_dawn_39","glia_P5_cake_16","P5_cake_40","glia_P5_duet_17","P5_duet_41","glia_P5_ebony_18","P5_ebony_42","glia_P6_alarm_19","P6_alarm_43","glia_P6_fog_20","P6_fog_44","glia_P6_gamma_21","P6_gamma_45","glia_P7_duck_22","P7_duck_46","glia_P7_east_23","P7_east_47","glia_P7_forest_24","P7_forest_48")
  # x = c("20200220_Mark_12836_1_1.1","20200220_Mark_12836_2_1.2","20200220_Mark_12836_3_1.3","20200220_Mark_12835_1_2.1","20200220_Mark_12835_2_2.2","20200220_Mark_12835_3_2.3","20200220_Mark_12835_4_2.4","20200220_Mark_12835_5_2.5","20200220_Mark_12570_1_3.1","20200220_Mark_12570_2_3.2","20200220_Mark_12570_3_3.3","20200220_Mark_12837_1_4.1","20200220_Mark_12837_2_4.2","20200220_Mark_12837_3_4.3","20200220_Mark_12837_4_4.4","20200220_Mark_12536_12_5.1","20200220_Mark_12536_16_5.2","20200220_Mark_12536_14_5.3","20200220_Mark_12536_15_5.4","20200220_Mark_12835_5_2.5replica")
  # x = c("KO 01","KO 02","KO 03", "WT 01","WT-02","WT-03","WT-04")
  # x = c("KO 02.1","KO 03.1","KO 04.2", "WT 01.2","WT-02","WT-03","WT-04")
  # cbind(x, metadata_matrix_from_filenames(x))
}



#' Count how often commonly used separation characters occur in an array of strings
#'
#' @importFrom stringr str_count
count_sep_char = function(x) {
  tibble(char_regex=c("_", " ", ";", "#", "-", "\\.", ",", "\\|")) %>%
    mutate(count = sapply(char_regex, function(s) sum(stringr::str_count(x, s))) ) %>%
    filter(count > 0) %>%
    arrange(desc(count))
}



# main split function
split_recursively = function(x, alt_sep_chars) {
  if(!is.matrix(x)) {
    return(x)
  }
  result = NULL
  for(j in 1:ncol(x)) {

    # possible improvement (niche cases, not high prio); regex replace just once in each element, starting regex replace from end of each string. if that results in a proper split, recurse
    j_split = split_char_num(x[,j])

    # if failed, let's try alternative sep chars
    if(!is.matrix(j_split)) {
      j_split = split_alt_sep_char(x[,j], sep_chars = alt_sep_chars)
    }

    # we split this column, so recurse to try and further subdivide
    if(is.matrix(j_split)) {
      j_split = split_recursively(j_split, alt_sep_chars)
    }

    result = cbind(result, j_split)
  }
  return(result)
}



# detect case of forgotten separation character. examples; wt1, wt2  ,  1plus, 1min, 2plus, 2min
split_char_num = function(x) {
  split_generic(x,
                list(c("([0-9])([a-z])", "\\1@\\2"),
                     c("([a-z])([0-9])", "\\1@\\2")),
                sep_char = "@")
}



# detect case of forgotten separation character. examples; wt.1, wt.2
split_alt_sep_char = function(x, sep_chars) {
  if(length(sep_chars) > 0) {
    return(split_generic(x, regex_list = list(c(paste0("([0-9a-z])(", paste(sep_chars, collapse="|"),")([0-9a-z])"), paste0("\\1@\\3") )), sep_char = "@"))
    # return(split_generic(x, list(c("([0-9a-z])[ .#-]([0-9a-z])", "\\1_\\2"))))
  }
  return(x)
}



score_classification_counts = function(x) {
  myclass = function(val) {
    cls = rep("char", length(val))
    cls[val == ""] = "na"
    cls[grepl("^[0-9]+$", val)] = "num"
    return(cls)
  }

  # higher table count = more repeated values = "better"
  # higher number of mixed data types (empty values ~ digits ~ characters) = not expected
  # myscore = function(x) mean(table(x)) / n_distinct(myclass(x)) # * mean(table(myclass(x)))
  # add punishment for inducing many empty values
  myscore = function(x) {
    x_cls = myclass(x)
    mean(table(x)) / n_distinct(x_cls)
  }

  if(is.matrix(x)) {
    return(sum(apply(x, 2, myscore)))
  } else {
    return(myscore(x))
  }

  ## default / simple case; only score counts
  # sum(apply(x, 2, function(x) mean(table(x))))
}



#' split_generic
#'
#' @importFrom stringr str_count str_split_fixed
split_generic = function(x, regex_list, sep_char = "@") {
  # insert split char between number and digit
  y = x
  for(r in regex_list) {
    y = gsub(r[1], r[2], y, ignore.case = T)
  }

  # if this transforms >90% of all elements in x...
  if(any(x != y)) { # if(sum(x == y) <= floor(length(x) * 0.1)) {
    # split and test if it yields repeated elements
    y_split = stringr::str_split_fixed(y, sep_char, n=Inf)

    # x_score = score_classification_counts(x)
    # y_score = score_classification_counts(y_split) / ncol(y_split)
    # if(x_score < y_score) {
    #   return(y_split)
    # }

    # count number of unique elements in each column -->> less unique elements = more certain there are legit repeated classification substrings
    x_count_ones = sum(table(x)==1)
    y_split_count_ones = apply(y_split, 2, function(x) sum(table(x)==1) )
    if(any(y_split_count_ones <= x_count_ones) ) {
      return(y_split)
    }
  }

  return(x)
}



# try to align matrix columns
# eg; inconsistent use of separation characters leads to misalignment. neuron_glia_1, something_1 -->> should have been neuron-glia_1
shift_empty = function(x) {
  if(!any(colSums(x=="") != 0)) return(x)
  # optionally, restrict to matrix with exactly 1 column that has empties
  # if(sum(colSums(x=="")) != 1) return(x)

  x_best = x
  x_best_score = score_classification_counts(x)
  ## shift from right to left
  for(j in ncol(x):2) {
    rows = x[,j] == ""
    if(any(rows)) {
      for(jj in (j-1):1) { # from current column j (containing empty values), shift to/from column jj
        # copy x, then shift into column j from the left (from column jj)
        x_new = x
        x_new[rows,(jj+1):j] = x_new[rows,jj:(j-1)]
        x_new[rows,jj] = ""
        x_score = score_classification_counts(x_new)
        if(x_best_score < x_score) {
          x_best = x_new
          x_best_score = x_score
        }
      }
    }
  }

  ## analogous, left to right
  for(j in 1:(ncol(x)-1)) {
    rows = x[,j] == ""
    if(any(rows)) {
      for(jj in (j+1):ncol(x)) { # from current column j (containing empty values), shift to/from column jj
        # copy x, then shift into column j from the right (up to column jj)
        x_new = x
        x_new[rows,j:(jj-1)] = x_new[rows,(j+1):jj]
        x_new[rows,jj] = ""
        # print(x_new)
        x_score = score_classification_counts(x_new)
        if(x_best_score < x_score) {
          x_best = x_new
          x_best_score = x_score
        }
      }
    }
  }

  return(x_best)
}



# Recognize dates formatted as 8 consecutive integers and rewrite them as "-" delimited
#
# If input array of strings are all encoding exactly 8 integers,
#   test if each represents a valid dates. If so, add separation chars.
#
# Implementation is intentionally simplified, only needs to catch most obvious cases (so ignore potential dates denoted as 6 numbers, etc.)
# Hardcoded "oldest year" is 1990, serves our use-case and adds robustness (we use this function in a scenario where any set of numbers can be supplied, but experiment samples are note likely 20+ years old)
# For ambiguous dates, e.g. 20102010 (20-10-2010 or 2010-20-10), prioritize year YYYYxxzz over xxzzYYYY
#
# examples;
# reformat_date(x = c("20200220", "20200221"))
# reformat_date(x = c("20202002", "20202102"))
# reformat_date(x = c("11012002", "11022002"))
reformat_date = function(x) {
  testYMD = function(y,md1,md2) {
    y = as.integer(y)
    md1 = as.integer(md1)
    md2 = as.integer(md2)
    all(y > 1990 & md1 >= 1 & md1 <= 12 & md2 >= 1 & md2 <= 31) ||
      all(y > 1990 & md1 >= 1 & md1 <= 31 & md2 >= 1 & md2 <= 12)
  }

  if(all(nchar(x) == 8) && all(grepl("^\\d+$", x))) {
    # YYYYMMDD or YYYYDDMM
    y = substr(x, 1, 4)
    md1 = substr(x, 5, 6)
    md2 = substr(x, 7, 8)
    if(testYMD(y, md1, md2)) {
      return(paste(y,md1,md2, sep="-")) # same order as input
    }
    # DDMMYYYY or MMDDYYYY
    y = substr(x, 5, 8)
    md1 = substr(x, 1, 2)
    md2 = substr(x, 3, 4)
    if(testYMD(y, md1, md2)) {
      return(paste(md1,md2, y, sep="-")) # same order as input
    }
  }
  return(x)
}
