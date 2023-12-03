from Bio import Entrez,SeqIO
import fasta_cross as fc


Entrez.email = 'david.bean@utahtech.edu'

# This gets multiple sequences from one organism
def get_sequences(organism_name,gene_name):

	Entrez.email = 'david.bean@utahtech.edu'

	search_query = f'{organism_name}[Organism] AND {gene_name}[Gene]'
	handle = Entrez.esearch(db ='protein',term = search_query)
	record = Entrez.read(handle)
	idlist = [i for i in record["IdList"]]

	records = {}
	for i in idlist:
		handle = Entrez.efetch(db = 'protein', id = i, rettype='gb',retmode='text')
		genbank_record = SeqIO.read(handle,'genbank')
		description = genbank_record.description+f' id:{i}'
		print(description+i)

		for feature in genbank_record.features:
			if feature.type == 'CDS':
				cds_feature = feature
				break
		sequence = cds_feature.location.extract(genbank_record.seq)
		print(sequence)

		records[description] = sequence
	return records





























database = 'protein'
def get_record(organism_name,gene_name):
	search_query = f'{organism_name}[Organism] AND {gene_name}[Gene]'
	handle = Entrez.esearch(db = database,term = search_query)
	record = Entrez.read(handle)

	return record

def get_genbank_record(i):
	handle = Entrez.efetch(db = database, id = i, rettype='gb',retmode='text')
	genbank_record = SeqIO.read(handle,'genbank')

	return genbank_record



def find_sequence(record):

	id_cds = None
	idlist = record['IdList']
	print(f'Len IdList = {len(idlist)}')

	for i in idlist:
		genbank_record = get_genbank_record(i)
		for feature in genbank_record.features:
			if feature.type == 'CDS':
				id_cds = i
				description = genbank_record.description+f' id:{i}'
				sequence = feature.location.extract(genbank_record.seq)
				break

	if id_cds == None:
		print(f'CDS not found in {genbank_record.description} genbank record')
	
	return [genbank_record.description,str(genbank_record)]




fc = fc.Fasta_cross()
fc.get_names(read_file = 'enterobacterales_short.fa',fasta_type = 'alex')
bacteria_names = fc.bacteria_names
gene_name = 'BamA'

for name in bacteria_names:
	record = get_record(organism_name = name, gene_name = gene_name)
	sequence = find_sequence(record)
	for thing in sequence:
		print(thing+'\n\n')











# with open('results.fasta','w') as f:
# 	for record in records:
# 		f.write(f'>{record[0]}\n{record[1]}\n')


# Only getting less than half back...

'''
	for i in idlist:
		if pull_sequence()

		for feature in genbank_record.features:
			if feature.type == 'CDS':
				cds_feature = feature
				break
'''