import fasta_filter as fasta_parsing
from sys import argv

if '-srt' in argv:
	sort_fasta = argv[argv.index('-srt')+1]
if '-nms' in argv:
	names_fasta = argv[argv.index('-nms')+1]

parser = fasta_parsing.Fasta_filter(sort_fasta,names_fasta)
parser.cross_ref_nrs()
parser.align()

