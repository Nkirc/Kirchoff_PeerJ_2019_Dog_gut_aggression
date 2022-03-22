# Kirchoff_PeerJ_2019_Dog_gut_aggression

This repository contains code used for the following publication:

Kirchoff NS, Udell MAR, Sharpton TJ. 2019. The gut microbiome correlates with conspecific aggression in a small population of rescued dogs (Canis familiaris) PeerJ 7:e6103 [https://doi.org/10.7717/peerj.6103] 

Data files are located on the Sharpton Lab Repository

Data files used in code located [here](http://files.cqls.oregonstate.edu/Sharpton_Lab/Papers/Kirchoff_PeerJ_2018/data/raw_16S/)

Raw 16S sequences located [here](http://files.cqls.oregonstate.edu/Sharpton_Lab/Papers/Kirchoff_PeerJ_2018/data/)

-----------------------------------------------------------------------

The code for this project is split into two R scripts. Claatu analysis and metacoder visualization code is in dog_claatu_analysis. All remaining code is found in Dog_Code_PeerJ

files needed for dog_claatu_analysis.R:

	- claatu_counts_rarefied.tab
	- dog_claatu_tax_dict.tab
	- mapping_file.txt
	- new_prepped_tree.tre
	- p_test_results.tab_stats.txt
	- rep_set_tax_assignments.txt

files needed for Dog_Code_PeerJ.R:

	-table_even40000.biom.txt
	-meta.3.csv
	-group_significance_Group.txt
	-unweighted_unifrac_test.df4_biom.txt
	-weighted_unifrac_test.df4_biom.txt

Raw 16S sequences are in the data directory. Metadata to match barcodes to samples are within files cg_primers_used_1_9_15_miseq and dog_1_15_map_file.txt.

#Dependencies

- R (v 3.2.3)
- vegan (v 2.3-3)
- coin (v 1.1-2)
- qvalue (v 2.2.2)

