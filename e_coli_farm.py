#! /usr/bin/env python3

import argparse
import re
import sys
from Bio import SeqIO

def farm_data(infile, target_species):
    records = list(SeqIO.parse(infile, "fasta"))
    target_species_index = None
    new_records = []
    for i, record in enumerate(records):
        if target_species in record.id:
            target_species_index = i
            break
    if target_species_index != None:
        target_sequence = records[target_species_index].seq
        # Create a new list of sequences with gaps removed for the target species

        for record in records:
            new_seq = ""
            for a, b in zip(target_sequence, record.seq):
                if a != "-":
                    new_seq += b
            new_record = record
            new_record.seq = new_seq
            new_records.append(new_record)
    return new_records

def generate_filename(infile):
    custom_filename = ""
    infile_without_ext = re.sub(r'\.(fasta|fa)$', '', infile)
    custom_filename = "{}_farmed.fa".format(infile_without_ext)
    return custom_filename

def write_farmed_fasta(outfile, records):
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
