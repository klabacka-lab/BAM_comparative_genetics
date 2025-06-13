from Bio import SeqIO
from Bio.Seq import Seq

def find_problematic_sequences(fasta_file):
    stop_codons = {"TAA", "TAG", "TGA"}
    problematic = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        stop_count = 0

        # Translate in-frame codons only
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in stop_codons:
                stop_count += 1

        if stop_count != 1:
            problematic.append(record.id)

    return problematic

# Example usage
if __name__ == "__main__":
    fasta_path = "working_bam.fa"
    results = find_problematic_sequences(fasta_path)
    print("Sequences with more than one or no stop codons:")
    for name in results:
        print(name)

