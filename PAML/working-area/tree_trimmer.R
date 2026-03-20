library(tidyverse)
library(ape)

## need to prune off a few branches, this is the file for 
tree = read.tree("PAML/working-area/unroot-tree.tree")
newtree = drop.tip(tree, tip=c("Hafnia_alvei","Rouxiella_chamberiensis","Erwinia_pyrifoliae")) 
print(newtree)
write.tree(newtree, file = "pruned_tree.tree")