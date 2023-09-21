# September 21 2023 (meeting notes with Alex)
# 1. Account for gaps for the proportion and unique values (gaps shouldn't be considered)
# 2. Create output file that provides the sites of highest conservation
# 3. Create sliding window of size 2 that moves one base at a time and takes the average of the proportion value
# 4. Once you do step 3 for sliding window of size 2, create a sliding window of varying sizes (you could have a function for sliding window that provides average proportion across a window size of n)
# 5. Identify two genes in Enterobacteriolis genomes. One that is predicted to be  highly conserved, and one that is predicted to be not highly conserved.  (maybe something David can do?)
# 6. Statistical tests to perform: An anova that compares the proportion of conservation across sites for the genes from #5 above, and an Anova that compares the conservation of sites between the different potra domains and the barrell itself.
# 7. Find out if Mycobacterium tuberculosis has a BamA homolog (take the M. tuberculosis genome, blast the e. coli BamA sequence against the genome, and see what gets hit) (maybe something David can do?)
# *** Alex- just start with steps 1-4 and work in order.

# September, 7 2023
# https://en.wikipedia.org/wiki/Terrabacteria
# Having a five part figure, where each part is a plot with conservation of BAMA as y axis and aa position as x axis
# Figures: see photos




# July 2023
# Method to read fasta file into dictionary using SeqIO:
## First, download biopython using pip or conda: https://biopython.org/wiki/Packages

from Bio import SeqIO
def read_fasta_file(file_path):
    sequence_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    return sequence_dict


