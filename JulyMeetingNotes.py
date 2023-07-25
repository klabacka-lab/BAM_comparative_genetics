# Method to read fasta file into dictionary using SeqIO:
## First, download biopython using pip or conda: https://biopython.org/wiki/Packages

from Bio import SeqIO
def read_fasta_file(file_path):
    sequence_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    return sequence_dict


