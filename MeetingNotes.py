# October 19 2023 (notes from meeting with Alex, David, and Randy)
## Create new alignment called 'bacteria.fa' that includes representatives of Fusobacteria, Aquificae, and Thermotogae
### Just a handful per lineage is fine
## Create new alignment called 'proteo_human.fa' 
### Include the human version of BAMA "
## Create new fastas for each of the alignments ('enterobacterales.fa', 'gammaproteo.fa', 'proteo.fa', 'hydrobacteria.fa', 'hydroterra.fa', 'bacteria.fa') where any position that doesn't have E. coli AA is removed
## To the csv files, add potra domains
## For next week, look into the BamA structure and map ultraconserved sites to the 3D structure.

# October 4 2023 (notes from meeting with David and Randy)
## It looks like SecY would be good as a conserved protein, and SecG would be good as a non-conserved protein. See following refs:
### SecY is conserved: https://www.sciencedirect.com/science/article/pii/S0005273602006624
### SecG does not show sequence conservation: https://www.nature.com/articles/367654a0
#### These two references were found here: https://journals.asm.org/doi/full/10.1128/jb.187.18.6454-6465.2005?casa_token=N3MrCJWvM4EAAAAA%3AIO1dUQTFCvXOy1d8lGEYCcOFmKNOhhbnfHN6eIphWcj9O68eD2oQMlfiAP17FMuxpRcllJ_0u4m1qYY



# October 2 2023 (notes from Randy)
# I spoke with Dr. Bakelar, and he recommended adding 'Helicobacter pylori' to our dataset
# He also recommended using SecY as our known conserved gene. As an alpha helical barrel complex that is also membrane-bound, it would be a good comparison for BamA.
# TamA might also be a good option, since it is similar to BamA (it also has potra domains)

# September 29 2923 (meeting notes from meeting with Alex and David)
# 1. Add a new column to the output.csv that is called "Domain". The values for each row should containt he domain that the position resides within (e.g., Potra1, Potra2, ... , BamA)

# September 21 2023 (meeting notes with Alex)
# 1. Account for gaps for the proportion and unique values (gaps shouldn't be considered)
# 2. Create output file that provides the sites of highest conservation
# 3. Create sliding window of size 2 that moves one base at a time and takes the average of the proportion value
# 4. Once you do step 3 for sliding window of size 2, create a sliding window of varying sizes (you could have a function for sliding window that provides average proportion across a window size of n)
# 5. Identify two genes in Enterobacteriolis genomes. One that is predicted to be  highly conserved, and one that is predicted to be not highly conserved.  (maybe something David can do?)
# 6. Statistical tests to perform: An anova that compares the proportion of conservation across sites for the genes from #5 above, and an Anova that compares the conservation of sites between the different potra domains and the barrell itself.
# 7. Find out if Mycobacterium tuberculosis has a BamA homolog (take the M. tuberculosis genome, blast the e. coli BamA sequence against the genome, and see what gets hit) (maybe something David can do?)
# *** Alex- just start with steps 1-4 and work in order.

# September, 7 2023
# https://en.wikipedia.org/wiki/Terrabacteria
# Having a five part figure, where each part is a plot with conservation of BAMA as y axis and aa position as x axis
# Figures: see photos




# July 2023
# Method to read fasta file into dictionary using SeqIO:
## First, download biopython using pip or conda: https://biopython.org/wiki/Packages

from Bio import SeqIO
def read_fasta_file(file_path):
    sequence_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    return sequence_dict


