#! /usr/bin/env python3

from Bio import Entrez
import argparse
import os
import sys
import time

def get_species(infile):
    with open(infile, "r") as infile:
        species_list = [line.strip() for line in infile.readlines()]
        return species_list

def search_genbank(gene, species, email):
    Entrez.email = email
    query = f"{gene} and {species}"
    handle = Entrez.esearch(db="nucleotide", term=query, idtype="acc", retmode="text")
    record = Entrez.read(handle)
    handle.close

    if record["IdList"]:
        return record["IdList"][0]
    else:
        return None

def retrieve_nuc_sequence(seq_ID, email):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=seq_ID, rettype="fasta", retmode="text")
    sequence = handle.read()
    handle.close()
    return sequence

def create_fasta_file(gene, species_list, tax_level, email):
    failure_count = 0
    outfile_path = f"{tax_level}_{gene}.fa"
    outfile = open(outfile_path, "w")
    failfile = open(f"{tax_level}_failed_species_{gene}.txt", "w")

    for species in species_list:
        ID_result = search_genbank(gene, species, email)

        if ID_result:
            sequence = retrieve_nuc_sequence(ID_result, email)
            outfile.write(sequence)
            cl_output = f"{species} sequence found and added | Total failures: {failure_count}"
            print(cl_output, end="\r")

        else:
            failure_count += 1
            failfile.write(f"{species} not found for {gene}" + "\n")
            cl_output = f"{species} sequence not found | Total failures: {failure_count}"
            print(cl_output, end="\r")

    outfile.close()
    failfile.close()
    print(f"{tax_level} FASTA file written.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
            type = str,
            help = "Input list of species of interest")
    parser.add_argument("-g", "--gene",
            type = str,
            help = "Gene of interest")
    parser.add_argument("-l", "--level",
            type = str,
            help = "Taxanomic level of interest (just used for naming files - no analysis)")
    parser.add_argument("-e", "--email",
            type = str,
            help = "Email address for use with Entrez to query NCBI")
    args = parser.parse_args()
    species_list = get_species(args.input)
    create_fasta_file(args.gene, species_list, args.level, args.email)

if __name__=="__main__":
    main()
