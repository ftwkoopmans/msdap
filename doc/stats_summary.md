
This vignette shows how to summarize MS-DAP DEA and/or differential
detection results in a dataset into a table with a single statistic per
gene and optionally map these to HGNC human gene identifiers.

## Summarizing statistical results from MS-DAP at gene level

After performing DEA, the result table from statistical analysis might
contain multiple proteingroups for a gene (i.e. isoforms). In downstream
analysis we often want to only condier only the most significant
proteingroup in such cases when performing GO analyses. Further, there
might be proteingroups that are ambiguous, i.e. multiple genes match the
uniprot accessions in the proteingroup (e.g. the gene_symbols_or_id
column shows “GRIA1;GRIA2”).

The `summarise_stats()` function helps filter/subset the DEA results
table to deal with this. Further, in case [differential
testing](differential_detection.md) has been applied to the dataset you
can also use this function to merge the respective top-hits into the
result table such that the output table contains 1 pvalue and effectsize
for each unique gene in your statistical results.

The `remove_ambiguous_proteingroups` parameter is critical, controlling
whether proteingroups that map to multiple genes, e.g. “GRIA1;GRIA2”,
are removed from the result table. This is an important decision; set to
`TRUE` to remove ambiguous proteingroups and return only proteins groups
with gene-unique evidence.

pro: this is the conservative option with lower odds of false positives
in downstream analyses, probably the better choice when using these
results to inform on followup experiments.

con: if there is no proteingroup with unique evidence for e.g. “GRIA1”
nor for “GRIA2”, removing the ambiguous “GRIA1;GRIA2” preludes
downstream GO analyses from finding enrichment for this gene (even
though for GO it wouldn’t have mattered whether we searched for GRIA1 or
GRIA2).

Note that in any case one shouldn’t search GO for *all* ambiguous genes
in a proteingroup (e.g. GRIA1 and GRIA2 in case of “GRIA1;GRIA2”)
without some form of correction because *either* protein was found, not
both. The ‘symbol’ column returned by this function will contain only
the ‘leading’ symbol; typically this column is used as input for GO (and
not the ‘gene_symbols_or_id’ which contains all/ambiguous symbols in
case `remove_ambiguous_proteingroups = FALSE`) So for GO analyses one
could retain ambiguous proteingroups and then follow the example below
to map these to unique human gene identifiers and retain 1 value per
unique gene

``` r
library(msdap)
# load your previously generated dataset from disk
load("C:/<path to your MS-DAP results>/dataset.RData")

# Summarise the statistical data in your dataset, assuming DEA was applied.
# note that if you analyzed multiple contrasts, optionally with multiple DEA algorithms, all results are appended in this table.
x = summarise_stats(
  dataset,
  # include DEA result in the result table
  return_dea = TRUE,
  # disabled in this example, but if differential detection was performed set this to TRUE to integrate respective 'strong z-scores'
  return_diffdetect = FALSE,
  # replace effectsizes with log2fc in result table. Niche usecase, e.g. optional for DEA models that apply shrinkage to foldchanges
  dea_logfc_as_effectsize = FALSE,
  # differential detection absolute z-score cutoff (only relevant when return_diffdetect=TRUE)
  diffdetect_zscore_threshold = 4,
  # should proteingroups that map to multiple genes, e.g. "GRIA1;GRIA2", be removed from the result table?
  remove_ambiguous_proteingroups = FALSE
)


# alternatively, when you've performed differential detection as well and want to merge the DEA results with the top-hits from differential detect
# first plot the differential detect histograms to double-check your desired z-score threshold
plot_differential_detect(dataset)
# compared to above example, we here set `return_diffdetect = TRUE`
x = summarise_stats(
  dataset, 
  return_dea = TRUE,                       
  return_diffdetect = TRUE,               
  dea_logfc_as_effectsize = FALSE,         
  diffdetect_zscore_threshold = 5,         # critical parameter
  remove_ambiguous_proteingroups = FALSE   # critical parameter
)
```

## Phylogenetic mapping from mouse proteins to human genes

The MS-DAP function `protein2gene_orthologs()` adds human orthologs for
mouse/rat proteins in a table that contains uniprot accessions and their
respective gene symbols.

This function is currently limited such that the first matching human
gene per proteingroup is returned, so ambiguous proteingroups with
multiple gene symbols (e.g. “GRIA1;GRIA2”) might not always yield the
‘leading’ gene and additional/ambiguous genes are disregarded). Columns
describing the HGNC (genenames.org) ID, NCBI Entez ID and Ensembl gene
ID are appended to the provided table.

Workflow for ID mapping, where subsequent steps are fallbacks;

1)  uniprot id -\> mgi/rgd id -\> hgnc id
2)  uniprot symbol -\> mgi/rgd id -\> hgnc id
3)  uniprot symbol -\> hgnc official symbol -\> hgnc id
4)  uniprot symbol -\> hgnc non-ambiguous synonyms -\> hgnc id

Before running below example code for phylogenetic mappings, download
required lookup tables here;

**HGNC**

- download link:
  <https://www.genenames.org/download/statistics-and-files/>
- table: “Complete dataset download links” –\>\> “Complete HGNC approved
  dataset text json” –\>\> download the “TXT” table
- filename is typically something like hgnc_complete_set.txt

**MGI**

- download link:
  <https://www.informatics.jax.org/downloads/reports/index.html>
- table: “MGI Marker associations to SWISS-PROT and TrEMBL protein IDs
  (tab-delimited)”
- filename is typically something like MRK_SwissProt_TrEMBL.rpt

for rat proteins one can use information from RGD as a drop-in
replacement for MGI in this example, see further `rgd_lookuptable()`.

``` r
# load HGNC and MGI lookup tables from disk, assuming you previously downloaded these
hgnc = hgnc_lookuptable(f = "C:/<path to your downloads>/hgnc_complete_set.txt")
mgi = mgi_lookuptable(f = "C:/<path to your downloads>/MRK_SwissProt_TrEMBL.rpt")

# map the mouse proteins to human gene identifiers
# note that the ID mapping success rate is printed to the console
x = protein2gene_orthologs(x, hgnc, mgi)

# print first few rows to console to inspect results
print(x, n=5)
# write full results to disk as Excel table
openxlsx::write.xlsx(x, "C:/<path to your MS-DAP results>/stats_summary.xlsx")

# If you retain ambiguous proteingroups in `summarise_stats()` and perform phylogenetic mapping, 
# some genes might have duplicate entries in your result table (e.g. "GRIA1" and "GRIA1;GRIA2" proteingroups will yield the same human gene ID).
# Here we'll remove all proteingroups that failed to map to human genes AND those that yield the same human ortholog gene.
# (this works because sorting of table x is already by 'best on top' after summarise_stats())
y = x %>% filter(!is.na(hgnc_id)) %>% distinct(hgnc_id, contrast, dea_algorithm, .keep_all = TRUE)
openxlsx::write.xlsx(y, "C:/<path to your MS-DAP results>/stats_summary_deduped.xlsx")
```

## Mapping from gene symbols directly to HGNC gene identifiers

For human proteome datasets we can directly match from uniprot gene
symbols to HGNC (genenames.org) identifiers by directly cross-matching
official symbols and synonyms. The function `protein2gene_by_symbol()`
can directly work with the ‘symbol’ column from the `summarise_stats()`
output. Do note that if you did not remove ambiguous proteingroups in
the `summarise_stats()` step, only the first gene symbol is used for ID
mapping here (for the subset of ambiguous proteingroups).

Before running below example code for phylogenetic mappings, download
the required lookup table here;

**HGNC**

- download link:
  <https://www.genenames.org/download/statistics-and-files/>
- table: “Complete dataset download links” –\>\> “Complete HGNC approved
  dataset text json” –\>\> download the “TXT” table
- filename is typically something like hgnc_complete_set.txt

``` r
# load HGNC lookup tables from disk, assuming you previously downloaded these
hgnc = hgnc_lookuptable(f = "C:/<path to your downloads>/hgnc_complete_set.txt")

# assuming x is the results you obtained from summarise_stats()
# this will add columns hgnc_id, hgnc_symbol, entrez_id and ensembl_id
x = protein2gene_by_symbol(x, hgnc)

# print first few rows to console to inspect results
print(x, n=5)

### below code is analogous to the mouse example above

# write full results to disk as Excel table
openxlsx::write.xlsx(x, "C:/<path to your MS-DAP results>/stats_summary.xlsx")

# If you retain ambiguous proteingroups in `summarise_stats()` and perform phylogenetic mapping, 
# some genes might have duplicate entries in your result table (e.g. "GRIA1" and "GRIA1;GRIA2" proteingroups will yield the same human gene ID).
# Here we'll remove all proteingroups that failed to map to human genes AND those that yield the same human ortholog gene.
# (this works because sorting of table x is already by 'best on top' after summarise_stats())
y = x %>% filter(!is.na(hgnc_id)) %>% distinct(hgnc_id, contrast, dea_algorithm, .keep_all = TRUE)
openxlsx::write.xlsx(y, "C:/<path to your MS-DAP results>/stats_summary_deduped.xlsx")
```
