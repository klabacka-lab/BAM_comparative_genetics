# We used the anvi'o platform to annotate bacterial genomes, extract genes of interest, create a multiple sequence alignment, and estimate a phylogenomic tree.
# Much of this script follows this tutorial: https://merenlab.org/2017/06/07/phylogenomics/
# Prodigal did give me some problems on my 2022 M2 mac, so I ran it on my 2020 M1 mac


# ------------------------------------------- #
# Step 0: set things up                       #
# ------------------------------------------- #

# First, install the anvi'o platform
## To install anvi'o, follow the instructions on their documentation: 
## https://anvio.org/install/
##  You may need to install prodigal separately; I followed their github documentation:
##  https://github.com/hyattpd/Prodigal/wiki/installation#installing-on-mac-os-x

# Rename all of the files so they end in 'fa'
for file in *.fa *.fna *.fasta; do

    # Extract the filename without its extension
    base_name=$(echo "$file" | sed -E 's/\.[^.]+$//')

    # Replace non-ASCII characters
	modified_base_name=$(echo "$original_var" | sed 's/[.() -]/_/g')

    # Create the new filename with the desired extension
    new_file="${modified_base_name}.fa"
	mv "$file" "$new_file"

	# Use anvi'o to reformat the files for processing
	anvi-script-reformat-fasta "$new_file" \
                           --simplify-names \
                           --overwrite-input
done


# ------------------------------------------- #
# Step 1: generate an anvi'o contigs database #
# ------------------------------------------- #

# create .db file for each .fa file
for i in `ls *fa | awk 'BEGIN{FS=".fa"}{print $1}'`
do
    anvi-gen-contigs-database -f $i.fa -o $i.db -T 4
    anvi-run-hmms -c $i.db
done

# ------------------------------------------- #
# Step 2: create contig database text file    #
# ------------------------------------------- #

# header row
echo -e "name\tcontigs_db_path" > external-genomes.txt

# add file base name and .db file name to database text file
ls *.db | awk -F'.db' '{print $1 "\t" $1 ".db"}' >> external-genomes.txt


# ------------------------------------------- #
# Step 3: see gene options for phylogenetics   #
# ------------------------------------------- #

# Look at all genes
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                   --hmm-source Bacteria_71 \
                                   --list-available-gene-names

# Look at 16S
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                   --hmm-source Ribosomal_RNA_16S \
                                   --list-available-gene-names

# ------------------------------------------- #
# Step 4: create concatenated alignment       #
# ------------------------------------------- #

# create concatenated alignment for all genes in the Bacteria_71 source
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                -o concatenated-genes.fa \
                                --hmm-source Bacteria_71 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate

# the new alignment 'concatenated-genes.fa' consists of the following genes:
#   ADK, AICARFT_IMPCHas, ATP-synt, ATP-synt_A,
#   Adenylsucc_synt, Chorismate_synt, EF_TS, Exonuc_VII_L, GrpE, Ham1p_like, IPPT,
#   OSCP, PGK, Pept_tRNA_hydro, RBFA, RNA_pol_L, RNA_pol_Rpb6, RRF, RecO_C,
#   Ribonuclease_P, Ribosom_S12_S23, Ribosomal_L1, Ribosomal_L13, Ribosomal_L14,
#   Ribosomal_L16, Ribosomal_L17, Ribosomal_L18p, Ribosomal_L19, Ribosomal_L2,
#   Ribosomal_L20, Ribosomal_L21p, Ribosomal_L22, Ribosomal_L23, Ribosomal_L27,
#   Ribosomal_L27A, Ribosomal_L28, Ribosomal_L29, Ribosomal_L3, Ribosomal_L32p,
#   Ribosomal_L35p, Ribosomal_L4, Ribosomal_L5, Ribosomal_L6, Ribosomal_L9_C,
#   Ribosomal_S10, Ribosomal_S11, Ribosomal_S13, Ribosomal_S15, Ribosomal_S16,
#   Ribosomal_S17, Ribosomal_S19, Ribosomal_S2, Ribosomal_S20p, Ribosomal_S3_C,
#   Ribosomal_S6, Ribosomal_S7, Ribosomal_S8, Ribosomal_S9, RsfS, RuvX, SecE,
#   SecG, SecY, SmpB, TsaE, UPF0054, YajC, eIF-1a, ribosomal_L24, tRNA-synt_1d,
#   tRNA_m1G_MT

# ------------------------------------------- #
# Step 5: extract 16S genes and blast         #
# ------------------------------------------- #

# BC01, BC03, BC04, and BC05 popped out at weird spots in the tree.
# I decided to extract 16S from each of these and blast it to find out species ID

anvi-get-sequences-for-hmm-hits -c BC01_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC01_draft_genome.16S.fa
# Results: no sequence obtained from this command; removed it from the dataset
anvi-get-sequences-for-hmm-hits -c BC03_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC03_draft_genome.16S.fa
# Results: Blasted top sequence in BC03_draft_genome.16S.fa, recovered Streptomyces globisporus at top.
## Total score: 16748, Query cover: 100%, E value: 0.0, % Ident: 100%
# Three of the top four hits were S. globisporus, with the only other being S. rubiginosohelvolus (at third)
## Total score: 16759, Query cover: 100%, E value: 0.0, % Ident: 100%
# We added the genomes of GCF_000718455_1_ASM71845v1_genomic (S. globisporus) and 
# GCF_014649875_1_ASM1464987v1_genomic (S. rubiginosohelvolus) to our dataset
anvi-get-sequences-for-hmm-hits -c BC04_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC04_draft_genome.16S.fa
# Results: Blasted top sequence in BC04_draft_genome.16S.fa, recovered "Acinetobacter radioresistens" at top.
## The top 8 species were all A. radioresistens, here are results for the top hit:
## Total score: 2824, Query cover: 100%, E value: 0.0, % Ident: 99.93%
# We added the genome GCF_003258335_1_ASM325833v1_genomic to our dataset    
anvi-get-sequences-for-hmm-hits -c BC05_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC05_draft_genome.16S.fa
# Results: Blasted top sequence in BC03_draft_genome.16S.fa, recovered "Streptomyces globisporus" as second.
## However, most of the top hits (including the top-most hit) are unnamed species (Streptomyces sp SM1P)
## S. globisporus is the second-ranked hit
## Total score: 2796, Query cover: 100%, E value: 0.0, % Ident: 100%
## All other described species have percent ident < 100%
## I then extracted the Ribosomal_RNA_23S and blasted it as well
## The top hit was also a Streptomyces ("S. sp") with 100% query cover, E value of 0.0, and percent identity of 99.55%.
## The second hit was Streptomyces parvus, with queyr cover of 100%, E value of 0.0, and percent identity of 99.58%.
## Even with added S. globisporus to the alignment, BC05 came out as sister to the Acinetobacter group (but highly divergent)
## I'm opting to remove it from the alignment as well- I think this may be due to contamination.

# Here are the blast hits for the other new genomes:
anvi-get-sequences-for-hmm-hits -c BC02_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC02_draft_genome.16S.fa
# Results: Blasted top sequence in BC02_draft_genome.16S.fa, recovered Streptomyces nigra at top.
## Total score: 16440, Query cover: 100%, E value: 0.0, % Ident: 99.54%
# Three of the top seven hits were S. nigra (one, two, and five), with the only other being S. sp
# However, there was not enough genomic data from this sample to include it in the phylogenomic analysis
## (The Bacteria_71 model in anvio had 0 hits)
anvi-get-sequences-for-hmm-hits -c BC06_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC06_draft_genome.16S.fa
# Results: Blasted top sequence in BC06_draft_genome.16S.fa, recovered "Streptomyces sp" at top.
## top 9 results are S. sp
## 10th result (with 99.27% identity, 100% query coverage) is S. collinus
## Added GCF_014204745_1_ASM1420474v1_genomic (S. collinus) to dataset
anvi-get-sequences-for-hmm-hits -c BC07_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC07_draft_genome.16S.fa
# Results: Blasted top sequence in BC07_draft_genome.16S.fa, recovered "Streptomyces luteogriseus" at top 2 spots.
## Total score: 16493, Query cover: 100%, E value: 0.0, % Ident: 99.6%
## Added GCF_014205055_1_ASM1420505v1_genomic to the dataset
anvi-get-sequences-for-hmm-hits -c BC08_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC08_draft_genome.16S.fa
# Results: Blasted top sequence in BC08_draft_genome.16S.fa, recovered "Streptomyces albogriseolus" at top spot.
# S. albogriseolus 4 of top 1 hits (all others S. sp)
## Total score: 2789, Query cover: 100%, E val: 0.0, % Ident: 99.93%
## NO REFSEQ GENOME
anvi-get-sequences-for-hmm-hits -c BC09_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC09_draft_genome.16S.fa
# Results: Blasted top sequence in BC09_draft_genome.16S.fa, recovered "Streptomyces nigra" at top spot.
## Total score: 16759, Query cover: 100%, E val: 0.0, % Ident: 100.00%
## Added GCF_003074055_1_ASM307405v1_genomic to dataset
anvi-get-sequences-for-hmm-hits -c BC10_draft_genome.db --hmm-source Ribosomal_RNA_16S -o BC10_draft_genome.16S.fa
# Results: Blasted top sequence in BC10_draft_genome.16S.fa, recovered "Streptomyces althioticus" and "S. lusitanus" at top spot.
## Total score: 16521, Query cover: 100%, E value: 0.0, % Ident: 100.00% (S. althioticus)
### NO REFSEQ GENOME FOR S. ALTHIOTICUS
## Total score: 16593, Query cover: 100%, E value: 0.0, % Ident: 100.00% (S. lusitanus)
## Added GCF_025427835_1_ASM2542783v1_genomic to dataset


# Additional taxa added to dataset for phylogenetics:
## GCF_000196555_1_ASM19655v1_genomic: Bifidobacterium longum
## GCF_011492945_1_ASM1149294v1_genomic: Rubrobacter tropicus
## GCF_001578075_1_ASM157807v1_genomic: Microcystis aeruginosa

# ------------------------------------------- #
# Step 6: estimate a phylogenomic tree        #
# ------------------------------------------- #

# Estimate a phylogeny from the concatenated dataset using the program "IQ-TREE"
## IQ-TREE multicore version 2.2.2.4
iqtree2 -s concatenated-genes.fa -m TEST -T AUTO -ntmax 18 -bb 1000