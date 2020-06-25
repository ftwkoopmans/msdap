#### public gist  @  https://gist.github.com/ftwkoopmans/7c4413811f08162ab35814b8259b09db
# code snippet to demonstrate MS-EmpiRe reproducibility issues & resolve by setting set.seed()
# msEmpiRe package build used for this testcase: https://github.com/zimmerlab/MS-EmpiRe/tree/ae985a92dada3ea6bce9ee88e6eb674ffe975f41
# author: Frank Koopmans, Department of Molecular and Cellular Neurobiology, VU University, Amsterdam, The Netherlands

library(msEmpiRe)

### example dataset and application of MS-EmpiRe following the documentation at: https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R
f <- system.file("extdata", "c1_c3.data", package = "msEmpiRe")
p <- system.file("extdata", "c1_c3.pdata", package = "msEmpiRe")
# loading data from installed data sets. The dataset is the yeast spike-in benchmarking data of O'Connell et al., as discussed in the paper
data <- msEmpiRe::read.standard(f, p,
                                prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                signal_pattern="c.*rep.*")
# extract the first two conditions
conditions <- msEmpiRe::extract_conditions(data)
conditions <- conditions[, c(1,2)]
# removing peptides that are detected in less than 2 samples per condition
eset_filtered <- msEmpiRe::filter_detection_rate(data, condition=conditions)
### end code snippet copied from documentation

## optional MS-EmpiRe normalization
# eset_filtered <- msEmpiRe::normalize(eset_filtered)


######### 1) run MS-EmpiRe twice, then compare results

# MS-EmpiRe differential expression analysis
msempire_result <- msEmpiRe::de.ana(eset_filtered)
print(table(msempire_result$p.adj <= 0.01)) # print the number of significant proteins

# on the issue of reproducibility; here we simply re-rerun MS-EmpiRe and check if results differ
msempire_result2 <- msEmpiRe::de.ana(eset_filtered)
print(table(msempire_result2$p.adj <= 0.01)) # print the number of significant proteins

# proteins are in the exact same order in result tables
# (thus both data tables are aligned and we can compare values on the same rows/indices)
stopifnot(msempire_result$prot.id == msempire_result2$prot.id)

# fold-changes are the same
fc_equals = all.equal(msempire_result$log2FC, msempire_result2$log2FC)
if(!isTRUE(fc_equals)) {
  warning("MS-EmpiRe analysis does not yield exact same foldchanges in a re-run of the same analysis")
  print(fc_equals)
  plot(msempire_result$log2FC, msempire_result2$log2FC)
}

# but the pvalue are not exactly the same
pval_equals = all.equal(msempire_result$p.val, msempire_result2$p.val)
if(!isTRUE(pval_equals)) {
  warning("MS-EmpiRe analysis does not yield exact same p-values in a re-run of the same analysis")
  print(pval_equals)
  plot(log10(msempire_result$p.adj), log10(msempire_result2$p.adj))
}



######### 2) fix reproducibility by hardcoding the random seed before each call to msEmpiRe::de.ana()

set.seed(123) # hardcoding the random seed
msempire_result <- msEmpiRe::de.ana(eset_filtered)
print(table(msempire_result$p.adj <= 0.01)) # print the number of significant proteins

set.seed(123) # hardcoding the random seed
msempire_result2 <- msEmpiRe::de.ana(eset_filtered)
print(table(msempire_result2$p.adj <= 0.01)) # print the number of significant proteins

# fold-changes and p-values are exactly the same
stopifnot(msempire_result$prot.id == msempire_result2$prot.id)
stopifnot(msempire_result$log2FC == msempire_result2$log2FC)
stopifnot(msempire_result$p.val == msempire_result2$p.val)



######### 3) confirm again that setting a distinct seed results in slightly altered p-values

set.seed(1234) # hardcoding the random seed
msempire_result3 <- msEmpiRe::de.ana(eset_filtered)
print(table(msempire_result3$p.adj <= 0.01)) # print the number of significant proteins

stopifnot(msempire_result$prot.id == msempire_result3$prot.id)
stopifnot(msempire_result$log2FC == msempire_result3$log2FC)
print(all.equal(msempire_result$p.val, msempire_result3$p.val))
stopifnot(msempire_result$p.val == msempire_result3$p.val)
