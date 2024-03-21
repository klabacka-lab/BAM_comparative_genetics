#! /usr/bin/env python3

from Bio import Entrez
import argparse
import os
import requests
import sys
import time

def get_species(infile):
    with open(infile, "r") as infile:
        species_list = [line.strip() for line in infile.readlines()]
        return species_list

def get_tax_ID(species, email):
    Entrez.email = email
    handle = Entrez.esearch(db="taxonomy", term=species)
    record = Entrez.read(handle)

    if record["Count"] != 0:
        return record["IdList"][0]
    else:
        return None

def species_to_tax_ID(species_list, email):
    tax_IDs = []
    start_time = time.time()
    progress_count = 0
    failure_count = 0

    for species in species_list:
        tax_ID = get_tax_ID(species, email)

        if tax_ID:
            tax_IDs.append(tax_ID)
            progress_count += 1
            current_time = time.time()
            elapsed_time = round(current_time - start_time)
            cl_output = f"Taxonomic IDs found: {progress_count} | Taxonomic IDs not found (error): {failure_count} | Time elapsed: {elapsed_time}"
            print(cl_output, end = "\r")

        else:
            failure_count += 1
            current_time = time.time()
            elapsed_time = round(current_time - start_time)
            cl_output = f"Taxonomic IDs found: {progress_count} | Taxonomic IDs not found (error): {failure_count} | Time elapsed: {elapsed_time}"
            print(cl_output, end="\r")

    return tax_IDs

def write_tax_ID_file(tax_IDs, outfile):
    with open(outfile, "w") as outfile:
        for tax_ID in tax_IDs:
            outfile.write(str(tax_ID) + "\n")
    print("List of taxonomic IDs written.")

def uniprot_search(gene, species):
    url = f"https://www.uniprot.org/uniprot/"
    parameters = {"query": f"gene:{gene}+organism_exact:\"{species}\"",
            "format": "fasta"}
    try:
        response = requests.get(url, params=parameters)
        if response.ok:
            fasta_data = response.text
            return fasta_data
        else:
            return None

    except requests.RequestException as e:
        print(f"Error: {e}")
        return None

def create_fasta_file(gene, species_list, tax_level):
    output_file = os.path.join("search_results", f"{gene}_{tax_level}.fasta")
    start_time = time.time()
    failure_count = 0

    for species in species_list:
        try:
            result = uniprot_search(gene, species)

            with open(output_file, "a") as outfile:
                outfile.write(result)

            current_time = time.time()
            elapsed_time = round(current_time - start_time)
            cl_output = f"{species} sequence found and added | Total failures: {failure_count} | Time elapsed: {elapsed_time}"
            print(cl_output, end="\r")

        except:
            failure_count += 1
            with open("failed_species.txt", "a") as out2file:
                out2file.write(f"{species} not found for {gene}" + "\n")

            current_time = time.time()
            elapsed_time = round(current_time - start_time)
            cl_output = f"{species} sequence not found | Total failures: {failure_count} | Time elapsed: {elapsed_time}"
            print(cl_output, end="\r")

    print(f"{tax_level} FASTA file written.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
            type = str,
            help = "Input list of species of interest")
    parser.add_argument("-t", "--taxa",
            type = str,
            default = "tax_IDs.txt",
            help = "Output list of taxanomic IDs")
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
#    tax_IDs = species_to_tax_ID(species_list, args.email)
#    print("\n")
#    write_tax_ID_file(tax_IDs, args.taxa)
#    print("\n")
    create_fasta_file(args.gene, species_list, args.level)

if __name__=="__main__":
    main()
