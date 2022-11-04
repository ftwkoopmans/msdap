
# Differential Expression Analysis

In this vignette we will apply MS-DAP to the LFQbench dataset and
compare the number of differentially expressed yeast proteins (spike-in,
foreground) to the human background proteins. Instead of applying the
entire MS-DAP pipeline which includes QC reporting etc., as seen in the
[user guide](userguide.md), this example demonstrates how to efficiently
use MS-DAP for a bioinformatic analysis by executing only functions
needed for DEA.

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

## classify yeast and human proteins

We here use a MS-DAP utility function to easily recognize the Human and
Yeast proteins, while ignoring E. Coli proteins and ambiguous
proteingroups (that match multiple of these classes). The
`regex_classification()` function applies the provided regular
expressions to the fasta headers of all proteins in each proteingroup to
determine classifications. If no fasta files were loaded prior to this,
as is the case in this example, the regular expressions are applied to
the protein identifiers instead.

``` r
dataset$proteins$classification = regex_classification(dataset$proteins$fasta_headers, regex=c(human="_HUMA", yeast="_YEAS", discard="_ECOL"))
print(table(dataset$proteins$classification))
#> 
#> discard   human   yeast 
#>    1312    3062    1786
```

## differential expression analyses

Below code demonstrates how to apply only the MS-DAP functions for
feature selection & normalization and DEA. Afterwards, we select only
human and yeast proteins from the DEA results in preparation of the next
analysis step. Note that the `analysis_quickstart()` function covers
most typical MS-DAP use-cases and provides lots of additional
functionality, this vignette shows how to specifically apply only
filtering and DEA by directly calling respective MS-DAP functions.

``` r
# feature selection: only peptides detected in 3+ replicates in each sample group, then apply normalization (vwmb algorithm, followed by between-group normalization at protein-level)
dataset = filter_dataset(dataset,
                         filter_min_detect = 3,
                         norm_algorithm = c("vwmb", "modebetween_protein"),
                         by_group = F, all_group = T, by_contrast = F)
#> progress: caching filter data took 2 seconds
#> progress: peptide to protein rollup with MaxLFQ (implementation: iq) took 1 seconds
#> info: filter dataset with settings: min_detect = 3; norm_algorithm = 'vwmb&modebetween_protein'; rollup_algorithm = 'maxlfq'
#> 12756/34263 peptides were retained after filtering over all groups
#> progress: peptide filtering and normalization took 2 seconds

# if you want to run the "msqrob" DEA algorithm instead of "ebayes" (and not use the analysis_quickstart() convenience function), you should first initialize multiprocessing by uncommending the following line;
# cl = initialize_multiprocessing(n_thread = 4)

# apply limma's eBayes to each contrast and flag proteins as significant at 5% FDR and foldchange larger than a threshold estimated from bootstrap analyses (specified by parameter; fc_signif=NA)
dataset = dea(dataset, dea_algorithm = "ebayes", qval_signif = 0.05, fc_signif = NA)
#> info: differential expression analysis for contrast: A vs B
#> info: using data from peptide filter: global data filter
#> progress: peptide to protein rollup with MaxLFQ (implementation: iq) took 1 seconds
#> info: log2 foldchange threshold estimated by bootstrap analysis: 0.652
#> progress: eBayes took 1 seconds

# add the yeast/human protein classifications to DEA score tibble and filter to only keep human and yeast proteins
tib_plot = left_join(dataset$de_proteins, dataset$proteins, by="protein_id") %>%
  filter(classification %in% c("human", "yeast"))
```

## print true/false-positives

Although the spike-in ratio in the LFQbench dataset is relatively large,
making this task easier than real-world datasets, this analysis
demonstrates that setting a cutoff for both the FDR corrected p-values
and the log2 foldchanges yields better results than only filtering by
p-value. The effect is larger for 5% FDR as compared to 1% FDR, which
naturally has fewer false-positives.

The foldchange threshold was estimated by a bootstrapping procedure
implemented in MS-DAP, refer to the *estimating foldchange thresholds*
section in the [introduction](intro.md) or source code of the
`dea_protein_background_foldchange_limits()` function for more details.

``` r
# 5% FDR  versus  5% FDR and foldchange threshold (taken together in column 'signif')
print(tib_plot %>% 
        group_by(classification, dea_algorithm) %>% 
        summarise(`5% FDR` = sum(qvalue <= 0.05),
                  `5% FDR AND foldchange threshold` = sum(signif)))
#> `summarise()` has grouped output by 'classification'. You can override using
#> the `.groups` argument.
#> # A tibble: 2 × 4
#> # Groups:   classification [2]
#>   classification dea_algorithm `5% FDR` `5% FDR AND foldchange threshold`
#>   <chr>          <chr>            <int>                             <int>
#> 1 human          ebayes             126                                16
#> 2 yeast          ebayes             629                               608

# analogous for more stringent q-value cutoff at 1% FDR
print(tib_plot %>% 
        group_by(classification, dea_algorithm) %>% 
        summarise(`1% FDR` = sum(qvalue <= 0.01),
                  `1% FDR AND foldchange threshold` = sum(qvalue <= 0.01 & signif)))
#> `summarise()` has grouped output by 'classification'. You can override using
#> the `.groups` argument.
#> # A tibble: 2 × 4
#> # Groups:   classification [2]
#>   classification dea_algorithm `1% FDR` `1% FDR AND foldchange threshold`
#>   <chr>          <chr>            <int>                             <int>
#> 1 human          ebayes              40                                10
#> 2 yeast          ebayes             599                               582
```
