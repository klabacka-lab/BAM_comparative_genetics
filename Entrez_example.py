from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO

Entrez.email = "randy.klabacka@utahtech.edu"

print("Hello! This script will download five sequences for your gene of interest")

gene_name = "HOXA1" #input("What is the gene name? ")

taxon1_name = "Amblyraja radiata"
taxon1_locus = "XM_033047630"

taxon2_name = "Rattus norvegicus" #input("What is the organism name for the second sequence? ")
taxon2_locus = "XM_017592997" #input("What is the locus name for the second sequence? ")

taxon3_name = "Cricetulus griseus"
taxon3_locus = "XM_003501983" #input("What is the locus name for the third sequence? ")

taxon4_name = "Ornithorhynchus anatinus" #input("What is the organism name for the fourth sequence? ")
taxon4_locus = "XM_029070115" #input("What is the locus name for the fourth sequence? ")

taxon5_name = "Xenopus tropicalis"
taxon5_locus = "NM_001008016" #input("What is the locus name for the third sequence? ")

taxon6_name = "Thamnophis elegans" #input("What is the organism name for the sixth sequence? ")
taxon6_locus = "XM_032234813" #input("What is the locus name for the sixth sequence? ")

taxon7_name = "Protopterus annectens" #input("What is the organism name for the sixth sequence? ")
taxon7_locus = "XM_044066791" #input("What is the locus name for the sixth sequence? ")

taxon8_name = "Homo sapiens" #input("What is the organism name for the first sequence? ")
taxon8_locus = "NM_005522" #input("What is the locus name for the first sequence? ")

taxon9_name = "Pan troglodytes"
taxon9_locus = "XM_016945561"

taxon10_name = "Mus musculus"
taxon10_locus = "NM_010449"

print("gene: " + gene_name)
print("organism 1: " + taxon1_name)
print("organism 2: " + taxon2_name)
print("organism 3: " + taxon3_name)
print("organism 4: " + taxon4_name)
print("organism 5: " + taxon5_name)
print("organism 6: " + taxon6_name)
print("organism 7: " + taxon7_name)
print("organism 8: " + taxon8_name)
print("organism 9: " + taxon9_name)
print("organism 10: " + taxon10_name)

def get_sequence(organism_name, gene_name, locus_name):
    search_query = f"{organism_name}[Organism] AND {gene_name}[Gene] AND {locus_name}[LOCUS]"
    handle = Entrez.esearch(db = "nucleotide", term=search_query)
    record = Entrez.read(handle)
    
    seq_id = record["IdList"][0]
    handle = Entrez.efetch(db = "nucleotide", id=seq_id, rettype="gb", retmode="text")
    genbank_record = SeqIO.read(handle, "genbank")
    for feature in genbank_record.features:
        if feature.type == "CDS":
            cds_feature = feature
            break
    cds_sequence = cds_feature.location.extract(genbank_record.seq)
    return cds_sequence

taxon1_seq = SeqRecord(get_sequence(taxon1_name, gene_name, taxon1_locus), id=taxon1_name)
taxon2_seq = SeqRecord(get_sequence(taxon2_name, gene_name, taxon2_locus), id=taxon2_name)
taxon3_seq = SeqRecord(get_sequence(taxon3_name, gene_name, taxon3_locus), id=taxon3_name)
taxon4_seq = SeqRecord(get_sequence(taxon4_name, gene_name, taxon4_locus), id=taxon4_name)
taxon5_seq = SeqRecord(get_sequence(taxon5_name, gene_name, taxon5_locus), id=taxon5_name)
taxon6_seq = SeqRecord(get_sequence(taxon6_name, gene_name, taxon6_locus), id=taxon6_name)
taxon7_seq = SeqRecord(get_sequence(taxon7_name, gene_name, taxon7_locus), id=taxon7_name)
taxon8_seq = SeqRecord(get_sequence(taxon8_name, gene_name, taxon8_locus), id=taxon8_name)
taxon9_seq = SeqRecord(get_sequence(taxon9_name, gene_name, taxon9_locus), id=taxon9_name)
taxon10_seq = SeqRecord(get_sequence(taxon10_name, gene_name, taxon10_locus), id=taxon10_name)

seq_records = [taxon1_seq, taxon2_seq, taxon3_seq, taxon4_seq, taxon5_seq, taxon6_seq, taxon7_seq, taxon8_seq, taxon9_seq, taxon10_seq]

SeqIO.write(seq_records, "unaligned.fasta", "fasta")

# Make sure you have mafft installed (see https://mafft.cbrc.jp/alignment/software/)
mafft_exe = "mafft"
mafft_cline = MafftCommandline(input="unaligned.fasta")
print(mafft_cline)
stdout, stderr = mafft_cline()
from Bio import AlignIO
with open("alignment.fasta", 'w') as handle:
    handle.write(stdout)

print(seq_records)
print("Alignment saved to 'alignment.fasta'")

# Estimate a phylogeny using maximum likelihood and bootstrapping through iqtree as a subprocess
import subprocess
input_alignment = "alignment.fasta"

# Define the IQ-TREE on for the subprocess
subprocess.run([
        "iqtree2",
        "-s", input_alignment,
        "-o", "Amblyraja",
        "-bb", "1000"
        ])

# Create image file from .treefile
# you may need to install ete3 from the command line, see https://anaconda.org/conda-forge/ete3

# import ete3 and relevant modules
from ete3 import Tree, TreeStyle

# load newick tree from file
tree_file = "alignment.fasta.treefile"
my_tree = Tree(tree_file)

# Render the tree to a png image
ts = TreeStyle()
ts.show_branch_support = True
output_png = "output_tree.png"
my_tree.render(output_png, w=800, h=600, dpi=300, tree_style=ts)




