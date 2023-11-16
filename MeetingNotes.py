# November 16 2023 (notes from meeting with Alex, David, and Randy)
## Today we got the analysis.r script running
## David is going to optimize the aligner using Entrez (currently running into some bugs with a few taxa)
## Alex will debug the analysis.r script (it is currently printing histograms with gaps)
## Alex will start drafting the poster (intro, methods, pipeline figure, playing with results slides)

# November 9 2023 (notes from meeting with Alex, David, and Randy)
## Alex created .csv files for all of the alignments ('enterobacterales.csv', 'gammaproteo.csv', 'proteo.csv', 'proteo_human.csv', 'hydrobacteria.csv', 'hydroterra.csv', 'bacteria.csv') for BamA. He hasn't done it for SecY because the formatting of the fasta wasn't conducive to the ecoli_farm.py script
## David is going to add a step to the crossreference/mining/aligning script (fastacross.py) that makes the naming convention conform to that needed for ecoli_farm.py
## David is going to look into the biopython Entrez package (to see if we can do the data mining directly from the command line instead of going through a web browser). Randy added a temporary file that uses Entrez (not eloquently) that can be used as a reference- it is called "Entrez_example.py"
## Alex is going to fuse the 'statistical_analysis.r' and 'conservation_histogram.r' files. The line that performs the analysis (lm) and the section that creates the box and whisker plot (everything under the "Filter out data with no domain" and "# Plot the data" comments) can be placed within the same function that creates the histograms.
## Alex will work on readme to integrate overview of bioinformatics pipeline (descriptions of each step, the scripts used, the file formats used and created), which David can convert into a illustrated version.
## IF YOU DARE: Look at figure 2a from the Copee paper we read this week- you can use it to inspire your adventorous soul

# November 2 2023 (notes from meeting with Alex, David, and Randy)
## Alex added domain information to the output files of stat.py
## We created the e_coli_farm.py script that removes sites from a fasta alignment that are not present in E. coli. We need to do this for SecY (conserved control) and another protein tbd (unconserved 'control')
## Randy made a draft r script to run an ANOVA that looks at the relationship between domain and conservation (with each position along the domain being a sample)
## The r script also plots the data (it is called 'statistical_analyses.r') and shows how you can create grids (joined plots in a single pdf)
## Alex is now going to create all of the csvs for BamA, SecY (after farmed), and unconserved_protein_tbd (after farmed) on all of the alignments ('enterobacterales.csv', 'gammaproteo.csv', 'proteo.csv', 'hydrobacteria.csv', 'hydroterra.csv', 'bacteria.csv') 


# October 26 2023 (notes from meeting with Alex, David, and Randy)
## See Oct 19 notes to create new fasta files for bacteria and proteo_human fasta
## Add domain names to csv outputs (using dictionary technique discussed in meeting)

# October 19 2023 (notes from meeting with Alex, David, and Randy)
## Create new alignment called 'bacteria.fa' that includes representatives of Fusobacteria, Aquificae, and Thermotogae
### Just a handful per lineage is fine
## Create new alignment called 'proteo_human.fa' 
### Include the human version of BamA (Sam50/Tob55)
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


