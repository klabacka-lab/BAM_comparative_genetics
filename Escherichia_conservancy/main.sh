#!/bin/bash

########################################
#### Retrieving Sequences###############
########################################
#
# IN: 	List of genes (ecoli_all_proteins)
#
#	Uniprot_search.py
#
# Out:	FASTA files, each containing sequences of one E.coli gene, and its homologous
#	counterparts in the other bacteria (if present)

# Uniprot_search.py assumes gene names will be in first column of tab deliniated file
#
# arg1:<gene file> arg2: <max number of genes to search>


geneFile="ecoli_all_proteins.txt"
maxNumGenes=5

python3 uniprot_search.py $geneFile $maxNumGenes 





#######################################
#### Filtering sequences###############
#######################################
#
# IN:	directory of FASTA files
#	array of species to keep
#	parameters to filter FASTA files:
#		min sequence lenght <= 50
#		contains sequences from all 5 species
#
#	sequence_filter.py
#
# OUT:	refined set of FASTA files


#######################################
#### Aligning sequences ###############
#######################################
#
# IN:	Unaligned FASTA files
#
#	aligner.py
#
# OUT:	Multiple sequence allignments


#######################################
#### farming alignments ###############
#######################################
#
# IN:	Multipe sequence alignments containing one gene from ecoli and its homologous
#	counterparts
#
#	e_coli_farm.py
#
# OUT:	Multiple sequence alignments in which insertaons not present have been trimmed


#######################################
#### Retrieving statistics ############
#######################################
#
# IN: 	cleaned multiple sequence alignments
#
#	position.py
#	stat.py
#
# OUT:	directory of CSV files (one per gene) containing proportion of most common
#	amino acid at each position of aligned sequences



#### Merging statistics
#
# IN: 	directory of CSV files
#	
#	Esch_stats.py (or R)
#
# OUT:	One CSV file containing all proteins and their average % conservancy
#	One CSV file containing all proteins and their median % conservancy



#### Placing BamA
#
# Output: of the genes in BamA, proportion that BamA is more conserved than


#### Histograms
#
#






#





#




#








#






#
