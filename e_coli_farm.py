#! /usr/bin/env python3

import argparse
import re
import sys
from Bio import SeqIO

def farm_gaps(infile, target_species):
    records = list(SeqIO.parse(infile, "fasta"))
    target_species_index = None
    for i, record in enumerate(records):
        if target_species in record.id:
            target_species_index = i
            break
    if target_species_index != None:
        target_sequence = records[target_species_index].seq
        new_records = [record for record in records if all(a == "-" or b != "-" for a,b in zip(target_sequence, record.seq))]
    return new_records

def generate_filename(infile):
    custom_filename = ""
    infile_without_ext = re.sub(r'\.(fasta|fa)$', '', infile)
    custom_filename = "{}_farmed.fa".format(infile_without_ext)
    return custom_filename

def write_farmed_fasta(outfile, records):
    with open(outfile, "w") as outhandle:
        SeqIO.write(records, outhandle, "fasta")
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
    farmed_records = farm_gaps(args.infile, args.species)
    write_farmed_fasta(outhandle, farmed_records)
