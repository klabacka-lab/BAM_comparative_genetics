from position import Position
from Bio import SeqIO
import csv
import re
import sys



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
        proportions.append(pos_obj.proportion)
    return proportions

def get_unique_aa(sites):
    unique_list = []
    for pos_obj in sites:
        unique_count = pos_obj.count_unique_aa()
        unique_list.append(unique_count)
    return unique_list

def get_sample_sizes(sites):
    sample_sizes_list = []
    for pos_obj in sites:
        sample_size = pos_obj.sample_size
        sample_sizes_list.append(sample_size)
    return sample_sizes_list

def excise_gaps(sites):
    for pos_obj in sites:
        sample_size = pos_obj.sample_size
        total = pos_obj.total_aa
        if sample_size / total < 0.07:
            pos_obj.proportion = "N/A"
    return sites

def find_conserved_sites(sites):
    conserved_sites = []
    for pos_obj in sites:
        if pos_obj != "N/A" and pos_obj.proportion > 0.98:
            site_aa = (pos_obj.pos, pos_obj.control_aa)
            conserved_sites.append(site_aa)
    return conserved_sites

def write_plot_csv(filename, positions, proportions, unique_counts, sample_sizes):
    custom_filename = generate_filename(filename, "stats")
    with open(custom_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Position', 'Proportion', 'Unique Amino Acids', 'Sample Sizes'])
        for pos, proportion, unique_count, sample_size in zip(positions, proportions, unique_counts, sample_sizes):
            csv_writer.writerow([pos + 1, proportion, unique_count, sample_size])

def write_conserved_csv(filename, sites):
    custom_filename = generate_filename(filename, "conservation")
    with open(custom_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Conserved Position', 'Proportion'])
        for pos in sites:
            csv_writer.writerow([pos])


def generate_filename(filename, purpose):
    filename_without_extension = re.sub(r'\.(fasta|fa)$', '', filename)
    custom_filename = ""
    if purpose == "stats":
        custom_filename = "{}_stats.csv".format(filename_without_extension)
    if purpose == "conservation":
        custom_filename = "{}_conserved.csv".format(filename_without_extension)
    return custom_filename

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(sys.argv[0] + ": Expecting alignment file path")
    fasta_file = str(sys.argv[1])
    seq_dict = read_fasta_file(fasta_file)
    sites_list_gaps = create_sites_list(seq_dict)
    sites_list = excise_gaps(sites_list_gaps)
    positions = list(range(len(list(seq_dict.values())[0])))
    proportions_list = get_proportions(sites_list)
    unique_counts_list = get_unique_aa(sites_list)
    sample_size_list = get_sample_sizes(sites_list)
    conserved_sites = find_conserved_sites(sites_list)
    write_plot_csv(fasta_file, positions, proportions_list, unique_counts_list, sample_size_list)
    write_conserved_csv(fasta_file, conserved_sites)


