#! /usr/bin/env python3

from Bio import Entrez, Medline

# Set your email (required by NCBI)
Entrez.email = "aeverett4746@gmail.com"

# Function to fetch accession numbers for a given species and gene
def fetch_accession_numbers(species_name, gene_name):
    # Construct the search term
    search_term = f"{gene_name}[Gene] AND {species_name}[Organism]"
    try:
        # Search the nucleotide database
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        # Return the list of accession numbers
        return record["IdList"]
    except:
        print(f"Oh, nice try, {species_name} wasn't found. Womp womp.")
# Read species names from a text file
with open("named_species.txt", "r") as file:
    species_list = [line.strip() for line in file.readlines()]

# Specify the gene name you're interested in
gene_name = "bamA"

# Initialize a dictionary to store results
accession_numbers = {}

# Fetch accession numbers for each species
for species in species_list:
    print(f"Searching for {gene_name} in {species}...")
    accession_numbers[species] = fetch_accession_numbers(species, gene_name)

# Print the results
with open("acc-nums.txt", "w") as outfile:
    for species, accessions in accession_numbers.items():
        if accessions:
            outfile.write(str(accessions[0]) + "\n")
        else:
            outfile.write(str(species))

print("Done-zo, Gonzo!")
