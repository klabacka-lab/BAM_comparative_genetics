'''
BamA protein database code: 6QGW
PyMOL
'''

from Bio.PDB import *
import os

pdbl = PDBList()
pdbl.retrieve_pdb_file('6QGW')

parser = MMCIFParser()
structure = parser.get_structure('6QGW','qg/6QGW.cif')




