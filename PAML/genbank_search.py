#! /usr/bin/env python3

from Bio import Entrez, SeqIO
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
    query = f"{species}[Organism] AND {gene}[Gene]"
    handle = Entrez.esearch(db="nucleotide", term=query, idtype="acc", retmode="text", retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        print(f"Found record for species: {species}")
        return record["IdList"][0]
    else:
        print(f"No record found for species: {species}")
        return None

def retrieve_gene_sequence(seq_ID, gene, email):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=seq_ID, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in record.features:
        if feature.type == "gene" and "gene" in feature.qualifiers and gene in feature.qualifiers["gene"]:
            gene_seq = record.seq[feature.location.start:feature.location.end]
            return f">{gene} {record.annotations['organism']} {seq_ID}\n{gene_seq}\n"
    
    print(f"Gene {gene} not found in {seq_ID}")
    return None

def create_fasta_file(gene, species_list, tax_level, email):
    failure_count = 0
    outfile_path = f"{tax_level}_{gene}.fa"
    outfile = open(outfile_path, "w")
    failfile = open(f"{tax_level}_failed_species_{gene}.txt", "w")

    for species in species_list:
        ID_result = search_genbank(gene, species, email)

        if ID_result:
            gene_sequence = retrieve_gene_sequence(ID_result, gene, email)
            if gene_sequence:
                outfile.write(gene_sequence)
                print(f"{species}: Gene sequence added | Total failures: {failure_count}", end="\r")
            else:
                failure_count += 1
                failfile.write(f"{species}: Gene not found in {ID_result}\n")
                print(f"{species}: Gene not found | Total failures: {failure_count}", end="\r")
        else:
            failure_count += 1
            failfile.write(f"{species}: No sequence found\n")
            print(f"{species}: No sequence found | Total failures: {failure_count}", end="\r")

    outfile.close()
    failfile.close()
    print(f"\n{tax_level} FASTA file written.")

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
            help = "Taxonomic level of interest (just used for naming files - no analysis)")
    parser.add_argument("-e", "--email",
            type = str,
            help = "Email address for use with Entrez to query NCBI")
    args = parser.parse_args()
    species_list = get_species(args.input)
    create_fasta_file(args.gene, species_list, args.level, args.email)

if __name__=="__main__":
    main()
