# BamA Comparative Genetics

The purpose of this project is to analyze the conservation of the protein BamA and its
satellite POTRA domains.

## How

It approaches this by looking in particular at the amino acid
sequences of BamA in Enterobacterales compared to other taxa of bacteria, as well as
against highly conserved and neutral genes in Enterobacterales bacteria.

## Methods

* position.py creates the Position class to calculate particular statistics at
given positions along the amino acid sequence alignment. Such statistics include
  * Sample Size : Number of amino acids present at given position, excluding
    gaps
  * Unique Amino Acids : The number of unique amino acids represented at given
    position
  * Proportion : The proportion of the most common amino acid to the total
    amount of amino acids found at given position

* stat.py takes a multiple sequence alignment FASTA file and tracks statistics
provided by the Position class across the entire alignment, writing them to CSV
files. Also removes significant gaps from alignments and makes note of positions
that are more highly conserved.


