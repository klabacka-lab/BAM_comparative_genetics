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


