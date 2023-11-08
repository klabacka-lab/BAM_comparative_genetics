#! /usr/bin/env python3

"""
A module that analyes the conservation of a multiple sequence amino acid alignment
"""

from position import Position
from Bio import SeqIO
import csv
import re
import sys

def read_fasta_file(file_path):

    """
    Reads a FASTA file and turns it into a manipulatable dictionary

    Paramaters
    ----------
    file_path: str
        File path for the FASTA file that is a multiple sequence amino acid alignment

    Returns
    -------
    sequence_dict: dict
        Dictionary composed of each sequence from FASTA file (Key = Sequence Descriptor, Value = Amino Acid Sequence)
    """

    sequence_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    return sequence_dict

def create_sites_list(seq_dict):

    """
    Creates list of position objects for each position on a multiple sequence alignment

    Parameters
    ----------
    seq_dict: dict
        Dictionary composed of each seqeuence from a FASTA file (Key = Sequence Descriptor, Value = Amino Acid Sequence)

    Returns
    -------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)
    """

    sites = []
    for pos_num in range(len(list(seq_dict.values())[0])):
        position = Position(pos_num, seq_dict)
        sites.append(position)
    return sites

def get_proportions(sites):

    """
    Collects proportion information from position objects
        Proportion = most common amino acid / total amino acids at a given position

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the postion on a multiple sequence alignment)

    Returns
    -------
    proportions: list
        List of proportion information for each position (each proportion being a float)
    """

    proportions = []
    for pos_obj in sites:
        proportions.append(pos_obj.proportion)
    return proportions

def get_unique_aa(sites):

    """
    Collects the number of unique amino acids across position objects

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)

    Returns
    -------
    unique_list: list
        List of unique amino acids at each position (each item being an integer representing total unique amino acids)
    """

    unique_list = []
    for pos_obj in sites:
        unique_count = pos_obj.count_unique_aa()
        unique_list.append(unique_count)
    return unique_list

def get_sample_sizes(sites):

    """
    Collects the sample size of each postion object
        Sample Size = total number of amino acids at a given position

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)

    Returns
    -------


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
        if pos_obj.proportion != "N/A" and pos_obj.proportion > 0.98:
            site_aa = (pos_obj.pos, pos_obj.control_aa)
            conserved_sites.append(site_aa)
    return conserved_sites

def write_plot_csv(filename, positions, proportions, unique_counts, sample_sizes, domain_labels):
    custom_filename = generate_filename(filename, "stats")
    with open(custom_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Position', 'Proportion', 'Unique Amino Acids', 'Sample Sizes', 'Domain'])
        for pos, proportion, unique_count, sample_size, domain_label in zip(positions, proportions, unique_counts, sample_sizes, domain_labels):
            csv_writer.writerow([pos + 1, proportion, unique_count, sample_size, domain_label])

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

def sliding_window(sites, window_size):
    averaged_proportions = []
    for i in range(len(sites) - window_size + 1):
        window = sites[i:i + window_size]
        total_prop = 0
        for position in window:
            if position.proportion != "N/A":
                total_prop += position.proportion
        average_prop = total_prop / window_size
        averaged_proportions.append(average_prop)
    return average_proportions

def label_domains(sites):
    domain_dict = {(24,91): "POTRA1", (92,172): "POTRA2", (175,263): "POTRA3", (266,344): "POTRA4", (347,421): "POTRA5", (448,810): "BamA"}
    domain_list = []
    for pos_obj in sites:
        labeled = False
        for start,stop in domain_dict.keys():
            if pos_obj.pos in range(start, stop+1):
                domain_list.append(domain_dict[(start,stop)])
                labeled = True
                break
        if not labeled:
            domain_list.append("")
    return domain_list

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
    domain_label_list = label_domains(sites_list)
    conserved_sites = find_conserved_sites(sites_list)
    write_plot_csv(fasta_file, positions, proportions_list, unique_counts_list, sample_size_list, domain_label_list)
    write_conserved_csv(fasta_file, conserved_sites)


