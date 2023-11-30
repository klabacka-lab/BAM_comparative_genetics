#! /usr/bin/env python3

"""
A module that "farms" an alignment for positions only present in a given species (default E.coli)
"""

import argparse
import re
import sys
from Bio import SeqIO

def farm_data(infile, target_species):

    """
    Farms alignment to remove gaps not aligned with target species

    Parameters
    ----------
    infile: str
        Filepath for multiple sequence alignment to be farmed

    target_species: str
        Species to target in the alignment

    Returns
    -------
    new_records: list
        List of SeqRecord objects with gaps farmed for target species
    """

    records = list(SeqIO.parse(infile, "fasta"))
    target_species_index = None
    new_records = []

    # Searches through FASTA file for target species
    for i, record in enumerate(records):
        if target_species in record.id and i != NONE:
            # Writes new farmed sequences based on target species
            target_sequence = records[i].seq
            new_seq = ""
            # Compares position in target sequence to current sequence
            for a, b in zip(target_sequence, record.seq):
                if a != "-":
                    new_seq += b
            new_record = record
            new_record.seq = new_seq
            new_records.append(new_record)
            break

    

def generate_filename(infile):

    """
    Generates custom output filename based on input filename

    Parameters
    ----------
    infile: str
        Filepath for input alignment

    Returns
    -------
    custom_filename: str
        Custom output filename for farmed alignment
    """

    custom_filename = ""
    infile_without_ext = re.sub(r'\.(fasta|fa)$', '', infile)
    custom_filename = "{}_farmed.fa".format(infile_without_ext)
    return custom_filename

def write_farmed_fasta(outfile, records):

    """
    Writes all farmed sequences to a new FASTA file

    Parameters
    ----------
    outfile: str
        Filepath for output FASTA file

    records: list
        List of farmed sequences to be written into new file

    Returns
    -------
    None
        Function writes new file
    """

    with open(outfile, "w") as outhandle:
        for record in records:
            outhandle.write(">" + record.id + "\n" + record.seq + "\n")
    print("Farmed file written.")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile",
            type = str,
            help = "FASTA file to be farmed.")
    parser.add_argument("-s", "--species",
            type = str,
            default = "Escherichiacoli",
            help = "Target species for farming.")
    args = parser.parse_args()
    outhandle = generate_filename(args.infile)
    farmed_records = farm_data(args.infile, args.species)
    write_farmed_fasta(outhandle, farmed_records)
