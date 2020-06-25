
  - [Bootstrapping procedure to estimate a threshold for
    foldchanges](#bootstrapping-procedure-to-estimate-a-threshold-for-foldchanges)
      - [load dataset](#load-dataset)
      - [classify yeast and human
        proteins](#classify-yeast-and-human-proteins)
      - [differential expression
        analyses](#differential-expression-analyses)
      - [print true/false-positives](#print-truefalse-positives)

# Bootstrapping procedure to estimate a threshold for foldchanges

To threshold Differential Expression Analysis (DEA) results by protein
foldchanges, users can optionally provide some cutoff value. But rather
than setting an arbitrary cutoff, you may also use MS-DAP to estimate a
threshold from N permutations of sample-to-condition assignments.
Permutations of sample labels within a group are disregarded as these
have no effect on the between-group foldchange, only unique combinations
of swapping samples between conditions A and B are considered. Finally,
we select the foldchange value at the 95% quantile of all permutations
as the threshold (`max(abs(quantile(fc_matrix, probs = c(1-probs,
probs), na.rm = T)))`). This is somewhat similar to the method described
by Hafemeister and Satija (<https://doi.org/10.1186/s13059-019-1874-1>).

## load dataset

1.  load the Skyline output of the LFQbench study (this file is bundled
    with the MS-DAP package as a test dataset).

2.  extract the respective group/condition of each sample by matching a
    regular expression against the filenames. This is just to
    demonstrate a more advanced convenience method that may come in
    handy for bioinformatic analyses, in typical MS-DAP workflows, the
    sample metadata is loaded from file (a csv or Excel table), see the
    ‘getting started’ vignette.

3.  finally, we define the contrast: group A versus group B.

note; for LFQbench dataset, the sample groups are taken from
process\_hye\_samples.R @
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
#> # A tibble: 6 x 2
#>   sample_id           group
#>   <chr>               <chr>
#> 1 lgillet_L150206_007 A    
#> 2 lgillet_L150206_009 A    
#> 3 lgillet_L150206_011 A    
#> 4 lgillet_L150206_008 B    
#> 5 lgillet_L150206_010 B    
#> 6 lgillet_L150206_012 B
```

## classify yeast and human proteins

Use MS-DAP utility function to easily recognize the human and yeast
proteins, while ignoring ecoli proteins or ambiguous proteingroups (that
match multiple of these classes).

``` r
dataset$proteins$classification = regex_classification(dataset$proteins$fasta_headers, regex=c(human="_HUMA", yeast="_YEAS", discard="_ECOL"))
print(table(dataset$proteins$classification))
#> 
#> discard   human   yeast 
#>    1312    3062    1786
```

## differential expression analyses

Note that these first two steps are incorporated in the quickstart
function, which should cover all common use-cases. We here demonstrate
how MS-DAP makes it easy to analyze any test dataset.

``` r
# feature selection: only peptides detected in 3+ replicates in each sample group, then apply normalization (vwmb algorithm)
dataset = filter_dataset(dataset,
                         filter_min_detect = 3,
                         norm_algorithm = "vwmb",
                         by_group = F, all_group = T, by_contrast = F)
#> progress: caching filter data took 2 seconds
#> progress: peptide filtering and normalization took 3 seconds

# apply limma's eBayes to each contrast and flag proteins as significant at 5% FDR and foldchange larger than a threshold estimated from bootstrap analyses (specified by parameter; fc_signif=NA)
dataset = dea(dataset, algo_de = "ebayes", qval_signif = 0.05, fc_signif = NA)
#> info: differential abundance analysis for contrast: A vs B
#> info: using data from peptide filter: global data filter
#> info: log2 foldchange threshold estimated by bootstrap analysis: 0.650
#> progress: eBayes took 1 seconds

# add the yeast/human protein classifications to DEA score tibble and filter to only keep human and yeast proteins
tib_plot = left_join(dataset$de_proteins, dataset$proteins, by="protein_id") %>%
  filter(classification %in% c("human", "yeast"))
```

## print true/false-positives

Although the spike-in ratio in the LFQbench dataset is relatively large,
making this task easier than real-world datasets, this example analysis
demonstrates that the additional foldchange thresholding mostly removes
false-positives from the list of significant hits. The effect is larger
for 5% FDR as compared to 1% FDR, which naturally has fewer
false-positives.

``` r
# 5% FDR  versus  5% FDR and foldchange threshold (taken together in column 'signif')
print(tib_plot %>% 
        group_by(classification) %>% 
        summarise(`5% FDR` = sum(qvalue <= 0.05),
                  `5% FDR AND foldchange threshold` = sum(signif)))
#> `summarise()` ungrouping output (override with `.groups` argument)
#> # A tibble: 2 x 3
#>   classification `5% FDR` `5% FDR AND foldchange threshold`
#>   <chr>             <int>                             <int>
#> 1 human               155                                16
#> 2 yeast               629                               603

# analogous for more stringent q-value cutoff at 1% FDR
print(tib_plot %>% 
        group_by(classification) %>% 
        summarise(`1% FDR` = sum(qvalue <= 0.01),
                  `1% FDR AND foldchange threshold` = sum(qvalue <= 0.01 & signif)))
#> `summarise()` ungrouping output (override with `.groups` argument)
#> # A tibble: 2 x 3
#>   classification `1% FDR` `1% FDR AND foldchange threshold`
#>   <chr>             <int>                             <int>
#> 1 human                48                                10
#> 2 yeast               596                               574
```
