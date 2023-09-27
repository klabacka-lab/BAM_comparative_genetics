from position import Position
from Bio import SeqIO
import csv
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
    print("Length of untampered sites: ", len(sites))
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

def excise_gaps(sites):
    for pos_obj in sites:
        if pos_obj.sample_size != pos_obj.total_aa:
            sites.remove(pos_obj)
    print("length of excised sites list: ", len(sites))
    return sites

def get_sample_sizes(sites):
    sample_sizes_list = []
    for pos_obj in sites:
        sample_sizes_list.append(pos_obj.sample_size)
    return sample_sizes_list

def get_total_aas(sites):
    total_list = []
    for pos_obj in sites:
        total_list.append(pos_obj.total_aa)
    return total_list

if __name__ == "__main__":
    fasta_file = "enterobacteralesOA.fa"
    seq_dict = read_fasta_file(fasta_file)
    sites_list = create_sites_list(seq_dict)
    gap_absent_sites = excise_gaps(sites_list)
    sample1 = get_sample_sizes(sites_list)
    sample1a = get_total_aas(sites_list)
    sample2 = get_sample_sizes(gap_absent_sites)
    sample2a = get_total_aas(gap_absent_sites)
    print("Number of sample sizes in untampered sites: ", len(sample1))
    print("Number of sample sizes in excised sites: ", len(sample2))
    print("Number of totals in untampered sites: ", len(sample1a))
    print("Number of totals in excised sites :", len(sample2a))
    print("Sample Size at Position 1 (including gaps?): ", sample1[0])
    print("Total Amino Acids at Postion 1 (including gaps?) ", sample1a[0])
    print("Sample Size at Position 1 (no gaps?) ", sample2[0])
    print("Total Amino Acids at Position 1 (no gaps?) ", sample2a[0])
