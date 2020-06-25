# devtools::load_all() # load our R package
#
# ## generate some test data
# npep = 4
# nsamp = 5
# nprot = 2
# samples = tibble(sample_id = paste0("s", 1:nsamp), group="grp1", exclude = F) %>% mutate(shortname=sample_id)
# peptides = tibble()
# for(j in sample(1:nprot, size = nprot)) { # randomize row order, to catch bugs that based on the peptide order
#   for(i in sample(1:npep, size = npep)) { # randomize row order, to catch bugs that based on the peptide order
#     # peptide fake log2 value; row number
#     # add some noise to each column, 0.1 extra to each additional column / sample / replicate
#     # selecting fixed value for this noise, we add more variation to the first (.1 relative to 1 is more noise than to subsequent rows), making the results of CoV based filtering nice and predictable. since values are in log2 space, multiply
#     x = rep(i, nsamp) * (1 + (0:(nsamp-1)) * 0.1)
#     #
#     peptides = bind_rows(peptides, tibble(protein_id=paste0("prot",j), peptide_id=paste0("prot",j,"_pep",i), sequence_plain=paste(rep(letters[i], 6), collapse=""), sequence_modified=paste(rep(letters[i], 6), collapse=""),
#                                           isdecoy=F, detect=T, rt=NA, sample_id=paste0("s", 1:nsamp), intensity=x))
#   }
# }
#
# # verify that after changing the detect counts, filtering works as expected
# # reducing the number of detects for a peptide decreases its score and subsequent filtering
# i_nondetect = which(grepl("prot2_pep3", peptides$peptide_id) & grepl("s[12]$", peptides$sample_id))
# peptides$detect[i_nondetect] = F
#
#
# # view the test data in wide format
# print( peptides %>% pivot_wider(id_cols = c(protein_id, peptide_id), names_from = "sample_id", values_from = "intensity") %>% arrange(protein_id, peptide_id) )
# print( peptides %>% pivot_wider(id_cols = c(protein_id, peptide_id), names_from = "sample_id", values_from = "detect") %>% arrange(protein_id, peptide_id) )
#
# # filter data
# DS = filter_dataset(cache_filtering_data(dataset = list(peptides=peptides, proteins=import_fasta(peptides), samples=samples)),
#                        filter_min_detect = 3, filter_fraction_detect = 0, filter_min_quant = 0, filter_fraction_quant = 0,
#                        filter_min_peptide_per_prot = 1,
#                        filter_topn_peptides = 3,
#                        norm_algorithm = "",
#                        by_group = T, all_group = T, by_contrast = F)
#
# # print pre-cached scores used in filtering, for debugging / QC
# DS$dt_pep_group %>% inner_join(DS$peptides %>% select(peptide_id, sample_id, key_peptide_group), by="key_peptide_group") %>% distinct(key_peptide_group, .keep_all = T) %>%
#   # drop 'key' columns and sort the table for easier reading/debugging
#   select(-starts_with('key_')) %>% arrange(peptide_id)
#
#
# ### print filtered data table
# # - without any difference in detect rates, in prot1, the worst CoV should be removed (pep4)
# # - introducing a difference in detect, which is reduced in 2 samples for prot2_pep3, the combined score (of detect&cov) now makes it worse than prot2_pep4
# #
# # - intensity_by_group  and  intensity_all_group  are the same in this case
# print( DS$peptides %>% pivot_wider(id_cols = c(protein_id, peptide_id), names_from = "sample_id", values_from = "intensity_by_group") %>% arrange(protein_id, peptide_id) )
# print( DS$peptides %>% pivot_wider(id_cols = c(protein_id, peptide_id), names_from = "sample_id", values_from = "intensity_all_group") %>% arrange(protein_id, peptide_id) )
#
