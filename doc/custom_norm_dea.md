
-   <a href="#load-dataset" id="toc-load-dataset">load dataset</a>
-   <a href="#custom-normalization-function"
    id="toc-custom-normalization-function">custom normalization function</a>
    -   <a href="#implementation" id="toc-implementation">implementation</a>
-   <a href="#custom-dea-function" id="toc-custom-dea-function">custom DEA
    function</a>
    -   <a href="#implementation-1" id="toc-implementation-1">implementation</a>
-   <a href="#lets-run-ms-dap-" id="toc-lets-run-ms-dap-">let’s run MS-DAP
    !</a>

This vignette demonstrates how you can plug custom functions into the
MS-DAP pipeline.

note: in the code blocks shown below (grey areas), the lines that start
with `#>` are the respective output that would be printed to the console
if you run these code snippets on your computer

## load dataset

1.  load the Skyline output of the LFQbench study (<PMID:27701404> \~
    this file is bundled with the MS-DAP package, you don’t have to
    download anything).

2.  extract the respective group/condition of each sample by matching a
    regular expression against the filenames. This is just to
    demonstrate a more advanced utility function that may come in handy
    for bioinformatic analyses. In typical MS-DAP workflows, the sample
    metadata would be loaded from file (a csv or Excel table), as
    demonstrated in the [user guide](userguide.md).

3.  finally, we define the contrast: group A versus group B.

note; the sample-to-condition assignments were taken from
`process_hye_samples.R` @
<https://www.ebi.ac.uk/pride/archive/projects/PXD002952/files>

``` r
library(msdap)

f <- system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")
dataset = import_dataset_skyline(f, confidence_threshold = 0.01, return_decoys = F, acquisition_mode = "dia")
#> info: reading Skyline report...
#> info: 4 unique target (plain)sequences ambiguously mapped to multiple proteins and thus removed. Examples; TTDVTGTIELPEGVEMVMPGDNIK, LNIISNLDCVNEVIGIR, LMDLSINK, EVDEQMLNVQNK
#> info: 34263/35943 precursors remain after selecting the 'best' precursor for each modified sequence

dataset = sample_metadata_custom(dataset, group_regex_array = c(A = "007|009|011", B = "008|010|012") )

dataset = setup_contrasts(dataset, contrast_list = list(c("A", "B")))
#> info: contrast: A vs B

print(dataset$samples %>% select(sample_id, group))
#> # A tibble: 6 × 2
#>   sample_id           group
#>   <chr>               <chr>
#> 1 lgillet_L150206_007 A    
#> 2 lgillet_L150206_009 A    
#> 3 lgillet_L150206_011 A    
#> 4 lgillet_L150206_008 B    
#> 5 lgillet_L150206_010 B    
#> 6 lgillet_L150206_012 B
```

## custom normalization function

We here implement a simple normalization approach to illustrate the
technique for providing MS-DAP with a custom normalization function. The
code below implements quantile normalization.

Note that in MS-DAP, normalization is always applied to peptide-level.
If you want to normalize at protein-level instead, we suggest to perform
peptide-to-protein rollup within the normalization function, apply
normalization and then backport those scaling factors to the
peptide-level matrix.

The requirements to make your custom normalization function work are:

1)  The function must have a unique name that does not overlap with
    common R functions (bad name: ‘mean’). So pick unique names like
    ‘mynorm_quantile’ or ‘my_norm_v1’

2)  It should at least have a parameter named `x_as_log2`. Furthermore,
    make sure to add a `...` parameter as a ‘catch all’ for named
    parameters passed by the MS-DAP pipeline that you won’t make use of,
    and to be robust against future parameters that might be added.

List of parameters that are currently passed to custom normalization
functions, add these to your function to work with the respective data:

-   `x_as_log2` = log2 transformed peptide intensity data matrix, where
    each row is a peptide and each column a sample
-   `group_by_cols` = array of characters that describe the sample
    grouping (which columns belong to the same sample group)
-   `group_by_rows` = array of characters that describe the peptide
    grouping (which rows belong to the same protein)
-   `rollup_algorithm` = user’s preferred peptide-to-protein rollup
    algorithm

3)  return data has to be of the *exact* same format as x_as_log2 (but
    with normalized abundance values)

### implementation

``` r
my_norm = function(x_as_log2, ...) {
  cat("Example custom normalization implementation, my_norm(), scaling all samples by some quantile\n")
  quantile_for_normalization = 0.95
  # value at quantile x for each column/sample
  scale_per_sample = matrixStats::colQuantiles(x_as_log2, probs = quantile_for_normalization, na.rm=T)
  # instead of scaling samples such that quantile x is zero, adjust by mean shift so output values are of the same order as input
  scale_per_sample = scale_per_sample - mean(scale_per_sample)
  # apply scaling to each column (transpose, scale 'rows', then transpose again)
  return(t(t(x_as_log2) - scale_per_sample))
}
```

now that we have implemented a normalization function, let’s run a very
basic data matrix to test if it works;

``` r
mock_data = cbind(1:6, 2:7, 0:5)
print(mock_data)
#>      [,1] [,2] [,3]
#> [1,]    1    2    0
#> [2,]    2    3    1
#> [3,]    3    4    2
#> [4,]    4    5    3
#> [5,]    5    6    4
#> [6,]    6    7    5
print(my_norm(mock_data))
#> Example custom normalization implementation, my_norm(), scaling all samples by some quantile
#>      [,1] [,2] [,3]
#> [1,]    1    1    1
#> [2,]    2    2    2
#> [3,]    3    3    3
#> [4,]    4    4    4
#> [5,]    5    5    5
#> [6,]    6    6    6
```

In the last section of this document, we’ll use this new normalization
function in MS-DAP

## custom DEA function

Below we have implemented a simplified DEA function that applies a
`t.test()` on protein abundances.

1)  The function must have a unique name that does not overlap with
    common R functions (bad name: ‘mean’). So pick unique names like
    ‘mydea_ttest’ or ‘my_ttest_v1’

2)  Analogous to the custom normalization function presented above,
    multiple named parameters will be passed to your custom dea
    function. You can choose which to work with, but make sure to add a
    `...` parameter as a ‘catch all’ for named parameters passed by the
    MS-DAP pipeline that you won’t make use of, and to be robust against
    future parameters that might be added.

List of parameters that are currently passed to custom dea functions,
add these to your function to work with the respective data:

-   `peptides` = long-format tibble of peptide data. For simple
    use-cases, you might first consider the peptide- or protein-level
    ExpressionSet (these hold simple data matrices to operate on, see
    example code block below)
-   `samples` = wide-format tibble with sample metadata, take note of
    the columns ‘sample_id’ and ‘condition’ so you can classify samples
    into conditions while A/B testing
-   `eset_peptides` = ExpressionSet for peptides. Holds a
    peptide\*sample log2 intensity matrix as well as row- and
    column-level sample data
-   `eset_proteins` = ExpressionSet for proteins, analogous to
    `eset_peptides` but peptide-to-protein rollup was already performed.
-   `input_intensities_are_log2` = boolean that tells you whether the
    intensities are log2 transformed (they are by default, but could
    change in future, so implement a check like below)
-   `random_variables` = array of column names in the samples table

The return data must be a tibble with these columns;

-   `protein_id` = character column with the protein identifiers
    provided in the input data (peptides tibble or metadata from
    peptide/protein ExpressionSet). These must be unique values
-   `pvalue` = numeric column with p-values
-   `qvalue` = numeric column with p-values adjusted for multiple
    testing correction (e.g. after applying FDR to your p-values)
-   `foldchange.log2` = numeric column with log2 foldchanges
-   `dea_algorithm` = name of your algorithm. must be a non-empty
    character string uniquely indicating the name/label of your method.
    DO NOT use names already used by other methods in this pipeline,
    like ‘ebayes’

### implementation

The example code, together with code comments, should provide a code
skeleton for the implementation of your custom DEA algorithm.

``` r
# note that we only add named parameters for data we're using in this example implementation  (i.e. additional parameters like `random_variables` and `peptides` are sent to `...`)
my_dea_stats = function(eset_proteins, input_intensities_are_log2, ...) {
  ### 1) from provided protein-level data; extract the protein intensity matrix, to which condition each sample belongs and find the columns matching groups 1 and 2
  x = Biobase::exprs(eset_proteins)
  # transform to log2 if input data is non-log
  if (!input_intensities_are_log2) {
    x = log2(x)
  }
  # set non-finite values to NA
  x[!is.finite(x)] = NA
  # extract groups
  x_groups = pData(eset_proteins)$condition
  x_groups_unique = unique(x_groups)
  # assume there are 2 groups
  stopifnot(length(x_groups_unique) == 2)
  cols_grp1 = x_groups == x_groups_unique[1]
  cols_grp2 = x_groups == x_groups_unique[2]

  ### 2) calculate p-values for each protein
  pval = rep(NA, nrow(x))
  for(i in 1:nrow(x)) {
    i_tt = t.test(x[i,cols_grp1], x[i,cols_grp2], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
    pval[i] = i_tt$p.value
  }
  # log2 fold-change
  # in MS-DAP, we compute B/A foldchanges (as noted in msdap::setup_contrasts() )
  log2FC = rowMeans(x[,cols_grp2,drop=F], na.rm=T) - rowMeans(x[,cols_grp1,drop=F], na.rm=T)

  ### 3) create a result tibble that contains all columns required for downstream compatability with this pipeline; protein_id, pvalue, qvalue, foldchange.log2, dea_algorithm
  result = tibble(protein_id = rownames(x),
                  pvalue = pval,
                  qvalue = p.adjust(pval, method = "fdr"),
                  foldchange.log2 = log2FC,
                  # the name of your algorithm. must be a non-empty character string uniquely indicating the name/label of your method
                  # do NOT use names already used by other methods in this pipeline, like 'ebayes' !
                  dea_algorithm = "my_dea_stats")

  cat("Example custom DEA implementation, my_dea_stats(), yielded", sum(is.finite(result$qvalue) & result$qvalue<=0.01), "hits at qvalue<=0.01\n")
  return(result)
}
```

## let’s run MS-DAP !

To use a custom function, we just write the *name of the R function* as
a parameter. Do not confuse this with the ‘dea_algorithm’ name given in
the result tibble of your dea function, that is merely the label used in
output plots and tables. We remember the LFQbench dataset is a
comparison of 2 conditions with 3 replicated measurements each, and
apply feature selection rules accordingly.

``` r
dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 2, # only peptides identified in at least 2 samples per sample group
  filter_min_quant = 3,  # only peptides quantified in at least 3 samples per sample group
  norm_algorithm = "my_norm", # custom algorithm !
  dea_algorithm = c("ebayes", "my_dea_stats"), # we use eBayes for reference, and also apply our custom dea algorithm independently !
  dea_qvalue_threshold = 0.01, # adjusted p-value to cutoff 'significant' hits @ dataset$de_proteins column 'signif'   (but note we don't use this output column in code below)
  dea_log2foldchange_threshold = 0, # foldchange threshold for proteins to be 'significant' @ dataset$de_proteins column 'signif' (but note we don't use this output column in code below)
  diffdetect_min_samples_observed = NA, # in this example, we disable the `differential detection` qualitative analysis
  output_dir = NA # in this example, we disable the generation of all output files
)
#> info: no output files will be generated, output_dir was set to NA
#> progress: caching filter data took 1 seconds
#> Example custom normalization implementation, my_norm(), scaling all samples by some quantile
#> Example custom normalization implementation, my_norm(), scaling all samples by some quantile
#> info: filter dataset with settings: min_detect = 2; min_quant = 3; norm_algorithm = 'my_norm'; rollup_algorithm = 'maxlfq'
#> 12756/34263 peptides were retained after filtering over all groups
#> 21591/34263 peptides were retained after filtering within each group independently ("by group")
#> progress: peptide filtering and normalization took 3 seconds
#> info: differential expression analysis for contrast: A vs B
#> info: using data from peptide filter: global data filter
#> progress: peptide to protein rollup with MaxLFQ (implementation: iq) took 1 seconds
#> progress: eBayes took 1 seconds
#> Example custom DEA implementation, my_dea_stats(), yielded 720 hits at qvalue<=0.01
```

Print the number of significant proteins

``` r
print(dataset$de_proteins %>% 
        group_by(dea_algorithm) %>% 
        summarise(`1% FDR` = sum(qvalue <= 0.01),
                  `1% FDR AND foldchange threshold` = sum(qvalue <= 0.01 & signif)))
#> # A tibble: 2 × 3
#>   dea_algorithm `1% FDR` `1% FDR AND foldchange threshold`
#>   <chr>            <int>                             <int>
#> 1 ebayes            1045                              1045
#> 2 my_dea_stats       720                               720
```

A simple Venn diagram of proteins at 1% FDR for each method shows our
custom function identifies a subset of the results from limma’s ebayes
function, so we verified that the example implementation doesn’t yield
random stuff.

``` r
gplots::venn(list(ebayes=dataset$de_proteins %>% filter(dea_algorithm=="ebayes" & qvalue <= 0.01) %>% pull(protein_id),
                  my_dea_stats=dataset$de_proteins %>% filter(dea_algorithm=="my_dea_stats" & qvalue <= 0.01) %>% pull(protein_id)))
```

![](images/dea-unnamed-chunk-8-1.png)<!-- -->
