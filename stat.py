from position import Position
from Bio import SeqIO
import csv

def read_fasta_file(file_path):
    sequence_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    return sequence_dict

def create_sites_list(seq_dict):
    sites = []
    for pos_num in range(len(list(seq_dict.values())[0])):
        position = Position(pos_num, seq_dict)
        sites.append(position)
    return sites

def get_proportions(sites):
    proportions = []
    for pos_obj in sites:
        prop = pos_obj.proportion()
        proportions.append(prop)
    return proportions

def get_unique_aa(sites):
    unique_list = []
    for pos_obj in sites:
        unique_count = pos_obj.unique()
        unique_list.append(unique_count)
    return unique_list

def get_sample_sizes(sites):
    sample_sizes_list = []
    for pos_obj in sites:
        sample_size = pos_obj.sample_size
        sample_sizes_list.append(sample_size)
    print(sample_sizes_list)
    return sample_sizes_list

def write_to_csv(filename, positions, proportions, unique_counts, sample_sizes):
    with open(filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Position', 'Proportion', 'Unique Amino Acids', 'Sample Sizes'])
        for pos, proportion, unique_count, sample_size in zip(positions, proportions, unique_counts, sample_sizes):
            csv_writer.writerow([pos + 1, proportion, unique_count, sample_size])

if __name__ == "__main__":
    fasta_file = "enterobacteralesOA.fa"
    seq_dict = read_fasta_file(fasta_file)
    sites_list = create_sites_list(seq_dict)
    positions = list(range(len(list(seq_dict.values())[0])))  
    proportions_list = get_proportions(sites_list)
    unique_counts_list = get_unique_aa(sites_list)
    sample_size_list = get_sample_sizes(sites_list)
    output_file = "entero_stats.csv"
    write_to_csv(output_file, positions, proportions_list, unique_counts_list)


