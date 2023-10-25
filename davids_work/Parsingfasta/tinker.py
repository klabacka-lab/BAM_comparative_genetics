import fasta_filter as ff
from sys import argv

if '-srt' in argv:
	sort_fasta = argv[argv.index('-srt')+1]
if '-nms' in argv:
	names_fasta = argv[argv.index('-nms')+1]
if '-out' in argv:
	write_handle = argv[argv.index('-out')+1]

parser = ff.Fasta_filter(sort_fasta,names_fasta,write_handle)
parser.cross_ref_nrs(write_handle = write_handle)
parser.align()

