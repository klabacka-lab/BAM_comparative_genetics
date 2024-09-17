from Bio import Entrez
from Bio import SeqIO

# Step 1: Set your email (required by NCBI)
Entrez.email = "d00369589@utahtech.edu"

# Step 2: Read species names from a .txt file
species_file = "named_species.txt"
with open(species_file, 'r') as f:
    species_names = [line.strip() for line in f]

no_genome_file = "no_genome_found.txt"
with open(no_genome_file, 'w') as ngf:
    ngf.write("Species with no genome found:\n")

# Step 3: Fetch genome data for each species
for species in species_names:
    print(f"Fetching genome for {species}...")

    # Step 4: Search GenBank for the species' genome
    search_term = f"{species}[SPECIES] AND complete genome[TITLE]"
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        # Step 5: Fetch the genome sequence using the first result
        genome_id = record["IdList"][0]
        fetch_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        try:
            genome_record = SeqIO.read(fetch_handle, "fasta")
            fetch_handle.close()

            # Step 6: Save the genome sequence to a GenBank file
            output_filename = f"{species.replace(' ', '_')}_genome.fasta"
            with open(output_filename, "w") as output_handle:
                SeqIO.write(genome_record, output_handle, "fasta")

            print(f"Saved genome for {species} as {output_filename}")
        except ValueError:
            print(f"No valid records found in GenBank format for {species}.")
            with open(no_genome_file, 'a') as ngf:
                ngf.write(f"{species}\n")
    else:
        print(f"No genome found for {species}")
        with open(no_genome_file, 'a') as ngf:
            ngf.write(f"{species}\n")

print("Process completed.")

