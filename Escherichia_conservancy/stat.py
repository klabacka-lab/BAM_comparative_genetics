#! /usr/bin/env python3

"""
A module that analyes the conservation of a multiple sequence amino acid alignment
"""

from position import Position
from Bio import SeqIO
import csv
import re
import sys
import os

def read_fasta_file(file_path):
    print('FILEPATH TO FASTA:',file_path)

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


    record = SeqIO.parse(file_path,"fasta")
    print(record)
    for record in SeqIO.parse(file_path, "fasta"):
        print(record)
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
        Sample Size = total number of amino acids at a given position (does not count gaps)

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)

    Returns
    -------
    sample_sizes_list: list
        List of sample sizes at each position (each item being an integer representing sample size)
    """

    sample_sizes_list = []
    for pos_obj in sites:
        sample_size = pos_obj.sample_size
        sample_sizes_list.append(sample_size)
    return sample_sizes_list

def excise_gaps(sites):

    """
    Negates positions that are more than 70% gap/missing amino acid

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)

    Returns
    -------
    sites: list
        Modified list of position objects with negated gaps
    """

    for pos_obj in sites:
        sample_size = pos_obj.sample_size
        total = pos_obj.total_aa
        if sample_size / total < 0.07:
            pos_obj.proportion = "N/A"
    return sites

def find_conserved_sites(sites):

    """
    Tracks position objects with a proportion value above 98% indicating high conservation

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about the position on a multiple sequence alignment)

    Returns
    -------
    conserved_sites: list
        List of position objects with high conservation of 98% or higher
    """

    conserved_sites = []
    for pos_obj in sites:
        if pos_obj.proportion != "N/A" and pos_obj.proportion > 0.98:
            site_aa = (pos_obj.pos, pos_obj.control_aa)
            conserved_sites.append(site_aa)
    return conserved_sites

def write_plot_csv(filename, positions, proportions, unique_counts, sample_sizes, domain_labels, directories = False):

    """
    Writes all collected data to a csv file for future plotting

    Parameters
    ----------
    filename: str
        Output filename for csv file

    positions: list
        List of integers indicating position on alignment

    proportions: list
        List of control proportions at each position (each proportion being a float)

    unique_counts: list
        List of unique amino acids at each position (each item being an integer)

    sample_sizes: list
        List of sample sizes at each position (each sample size being an integer)

    domain_labels: list
        List of domain labels (strings indicating which domain each position is located under)

    Returns
    -------
    None
        Function writes a csv file
    """

    custom_filename = generate_filename(filename, "stats",directories = False)
    with open(custom_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Position', 'Proportion', 'Unique Amino Acids', 'Sample Sizes', 'Domain'])
        for pos, proportion, unique_count, sample_size, domain_label in zip(positions, proportions, unique_counts, sample_sizes, domain_labels):
            csv_writer.writerow([pos + 1, proportion, unique_count, sample_size, domain_label])
    print("Stat file created.")

def write_conserved_csv(filename, sites, directories):

    """
    Writes most conserved sites of the alignment to a csv file for future plotting

    Parameters
    ----------
    filename: str
        Output filename for csv file

    sites: list
        List of position objects for sites conserved at 98% or above

    Returns
    -------
    None
        Function writes a csv file
    """

    custom_filename = generate_filename(filename, "conservation", directories)

    with open(custom_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Conserved Position', 'Amino Acid'])
        for pos in sites:
            csv_writer.writerow([pos])
    print("Conserved site file written.")

def generate_filename(filename, purpose, directories):

    """
    Generates a custom filename for csv files based on functionality

    Parameters
    ----------
    filename: str
        Input filename to be manipulated

    purpose: str
        Determines function of csv file for custom filename

    Returns
    -------
    custom_filename: str
        Custom output filename specified for function
    """

    if directories == False:

        filename_without_extension = re.sub(r'\.(fasta|fa)$', '', filename)
        custom_filename = ""
        if purpose == "stats":
            custom_filename = "{}_stats.csv".format(filename_without_extension)
        if purpose == "conservation":
            custom_filename = "{}_conserved.csv".format(filename_without_extension)
        return custom_filename

    else:

        filename_without_extension = re.sub(r'\.(fasta|fa)$', '', filename)
        custom_filename = ""
        if purpose == "stats":
            custom_filename = "./stats/"+"{}_stats.csv".format(filename_without_extension)
        if purpose == "conservation":
            custom_filename = "./conserved/"+"{}_conserved.csv".format(filename_without_extension)
        return custom_filename

def sliding_window(sites, window_size):

    """
    Come back to this, Alex
    """

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

    """
    Labels position objects by domain name based on position location on alignment

    Parameters
    ----------
    sites: list
        List of position objects (each object containing information about a position on the multiple sequence alignment)

    Returns
    -------
    domain_list: list
        List of labels for each position based on location in reference to domains (based on Escherichia coli annotated protein)
    """

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


    fasta_file = str(sys.argv[1]) #command line input filepath
    print(os.getcwd())



    # David's work to get it working within Esch pipeline
    # Changing working directoryt to farmed_results
    fasta_file = fasta_file.split('/')[1]
    os.chdir('farmed_results')
    print(os.getcwd())



    seq_dict = read_fasta_file(fasta_file)
    print(seq_dict)
    sites_list_gaps = create_sites_list(seq_dict) #sites with no alteration
    sites_list = excise_gaps(sites_list_gaps) #updated sites list with significant gaps labeled
    positions = list(range(len(list(seq_dict.values())[0]))) #generates list of position numbers
    proportions_list = get_proportions(sites_list)
    unique_counts_list = get_unique_aa(sites_list)
    sample_size_list = get_sample_sizes(sites_list)
    domain_label_list = label_domains(sites_list)
    conserved_sites = find_conserved_sites(sites_list)



    # Changing working directory from farmed_results to stats
    os.chdir('../stats')
    print(os.getcwd())
    write_plot_csv(fasta_file, positions, proportions_list, unique_counts_list, sample_size_list, domain_label_list,directories = True)

    # Changing working directory from stats to conserved
    os.chdir('../conserved')
    print(os.getcwd())
    write_conserved_csv(fasta_file, conserved_sites,directories = False)


