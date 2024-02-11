#!/bin/bash

########################################
#### Retrieving Sequences###############
########################################
#
# IN: 	List of genes (ecoli_all_proteins)
#	taxID to limit search
#
#	Uniprot_search.py
#
# Out:	FASTA files, each containing sequences of one E.coli gene, and its homologous
#	counterparts in the other bacteria (if present)

# Uniprot_search.py assumes gene names will be in first column of tab deliniated file
#
# arg1:<gene file> arg2: <max number of genes to search>


geneFile="ecoli_all_proteins.txt"
maxNumGenes=10
taxID=561
echo ""
echo "==================================================================="
echo "== STEP 1: Retrieving Sequences ==================================="
echo "==================================================================="
echo ""
echo "Running uniprot searches for genes in $geneFile within taxID $taxID"
echo "Limiting search to the first $maxNumGenes genes in $geneFile"
echo ""
# clearing previous search results

if [ -d "search_results" ]; then
	echo "Removing previous search results"
	rm -r search_results
fi

if [ -d "filtered_results" ]; then
	rm -r filtered_results
fi

if [ -d "farmed_results" ]; then
	rm -r farmed_results
fi 


mkdir search_results
python3 uniprot_search.py $geneFile $taxID $maxNumGenes 
echo ""

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
# arg1 <minimum sequence length> arg2 <array of species of interest>

# Note: for speciesOfInterest array, delimit species with a comma. Use underspaces instead
# 	of spaces.
# Potential future change: remove search results directory to save space/memory

minSeqLen=50
speciesOfInterest=("Escherichia_albertii,Escherichia_coli,Escherichia_fergusonii,Escherichia_marmotae,Escherichia_ruysiae")

mkdir filtered_results
echo "==================================================================="
echo "== STEP 2: Filtering Results ======================================"
echo "==================================================================="
echo ""
echo "Filtering retrieved sequences to following specifications:"
echo "	minimum sequence length: $minSeqLen"
echo "	must contain: ${speciesOfInterest[@]}"
echo ""
python3 sequence_filter.py $minSeqLen ${speciesOfInterest[@]}
echo ""



#######################################
#### Aligning sequences ###############
#######################################
#
# IN:	Unaligned FASTA files
#
#	aligner.py
#
# OUT:	Multiple sequence allignments
# NOTE: muscle must be installed for this step
echo "==================================================================="
echo "== STEP 3: Aligning ==============================================="
echo "==================================================================="
echo ""
echo "NOTE: for this step, muscle must be installed."
echo ""
python3 aligner.py
echo ""
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
echo "==================================================================="
echo "== STEP 4: Farming  ==============================================="
echo "==================================================================="
mkdir farmed_results
echo ""
infileDir=filtered_results
outfileDir=farmed_results
farmSpecies="Escherechia coli"


for filename in $infileDir/*.fasta; do
	#python3 e_coli_farm.py -i $filename -s $farmSpecies -d $outfileDir
	python3 e_coli_farm.py -i $filename -d $outfileDir -q true
done
echo "Farming complete"
echo""

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
