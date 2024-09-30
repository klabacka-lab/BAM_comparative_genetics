import csv
from Bio import SeqIO

genomes_dir = './collection/genomes/fasta/'
output_file = 'bacteria_bamA.fa'

species_positions = []
with open('species_extract.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        species_positions.append({
            'species': row['species'],
            'start': int(row['start']),
            'stop': int(row['stop'])
        })

with open(output_file, 'a') as out_fasta:
    for entry in species_positions:
        species = entry['species']
        start = entry['start']
        stop = entry['stop']
        genome_file = f'{genomes_dir}{species}_genomic.fa'
        with open(genome_file, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                extracted_seq = record.seq[start-1:stop]
                out_fasta.write(f'>{species}_{start}_{stop}\n{extracted_seq}\n')
