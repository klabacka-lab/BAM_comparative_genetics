# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import ape library
## Install if you haven't using install.packages('ape'=Analysis of Phylogenetics and Evolution)
library(ape)
library(phytools)

# Read tree (which is a Newick file)
bact_tree <- read.tree('opentree_relationships.tre')

# Prune branches that were not placed
## I put "REMOVE" for tips without scientific name using regex in vim
### :%s/mrca\w*/REMOVE/g
## I removed ott numbers in vim
### :%s/\([a-zA-Z]*\)_ott\d*/\1/g
tips_to_remove <- c("REMOVE","Escherichia_marmotae", "Shimwellia_blattae", "Serratia_entomophila", "Biostraticola_tofi", "Plesiomonas_shigelloides", "Cedecea_davisae", "Buttiauxella_ferragutiae", "Buttiauxella_brennerae", "Buttiauxella_warmboldiae", "Buttiauxella_izardii", "Buttiauxella_noackiae", "Buttiauxella_agrestis", "Morganella_morganii", "Raoultella_electrica", "Kluyvera_ascorbata", "Orientia_tsutsugamushi")
bact_tree_pruned <- drop.tip(bact_tree, tips_to_remove)

# get tip labels
tips <- bact_tree_pruned$tip.label
writeLines(tips, "tree_tips.csv")

bact_tree_pruned <- ladderize(bact_tree_pruned)

plotTree(bact_tree_pruned, type="fan", fsize = 0.6, edge.width = 2, label.offset = 0.5)



