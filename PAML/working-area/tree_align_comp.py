''' planning
load tree file
read name from tree file individually


'''
import re


def prep_tree(tree_file):
    with open(tree_file, "r") as tree:
        newtree = []
        tree_txt = tree.read().split(",")
    for element in tree_txt:
        element = element.split(':')[0].lstrip("(((((").split("_")
        newtree.append("_".join(element[0:2]))
    return newtree

def prep_align(align_file):
    with open(align_file, "r") as align:
        align_list = []
        align_txt = align.read().split("\n")
        for element in align_txt:
            align_list.append(element.split(" ")[0])
        return align_list

def compare(new_tree, align_list):
    retlist = []
    for element in new_tree:
        if element not in align_list:
            retlist.append("tree " + element)
    for element in align_list:
        if element not in new_tree:
            retlist.append("align " + element)
    return retlist


######


print(compare(prep_tree(input("Tree: ")), prep_align(input("Align: "))))