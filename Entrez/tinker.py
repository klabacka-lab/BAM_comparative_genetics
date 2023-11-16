from Bio import Entrez,SeqIO
import fasta_cross as fc
Entrez.email = "david.bean@utahtech.edu"


'''
cds_feature error: Buttiauxella ferragutiae
cds_feature error: Buttiauxella izardii
cds_feature error: Buttiauxella noackiae
cds_feature error: Buttiauxella warmboldiae
cds_feature error: Cedecea davisae
cds_feature error: Citrobacter gillenii
cds_feature error: Cronobacter malonaticus
cds_feature error: Cronobacter muytjensii
cds_feature error: Cronobacter sakazakii
cds_feature error: Cronobacter turicensis
cds_feature error: Enterobacter cloacae
cds_feature error: Escherichia hermannii
cds_feature error: Escherichia marmotae
cds_feature error: Escherichia ruysiae
cds_feature error: Klebsiella aerogenes
cds_feature error: Klebsiella huaxiensis
cds_feature error: Klebsiella michiganensis
record['IdList'][0] out of range error:Phytobacter palmae

tried bamA. 42 sequences returned. Trying BamA. Same. capitalization does not seem to matter
'''

organism_name = 'Buttiauxella ferragutiae'
gene_name = 'bamA'

query = f'{organism_name}[Organism] AND {gene_name}[GENE]'
handle = Entrez.esearch(db='protein',term= query)
record = Entrez.read(handle)
IdList = record["IdList"]

handle = Entrez.efetch(db = 'protein', id = IdList[1], rettype='gb',retmode='text')
genbank_record = SeqIO.read(handle,"genbank")

for feature in genbank_record.features:
	print(feature)