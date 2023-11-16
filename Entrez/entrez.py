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


# Gets the first sequence and description for a given organism
def get_sequence(organism_name,gene_name):
	search_query = f'{organism_name}[Organism] AND {gene_name}[Gene]'
	handle = Entrez.esearch(db ='protein',term = search_query)
	record = Entrez.read(handle)



	if len(record["IdList"]) != 0:
		ident = record["IdList"][0]
	else:
		print(f"record['IdList'][0] out of range error:{organism_name}")
		return None

	handle = Entrez.efetch(db = 'protein', id = ident, rettype='gb',retmode='text')
	genbank_record = SeqIO.read(handle,'genbank')
	description = genbank_record.description+f' id:{ident}'

	index = 0
	cds_feature = None

# id from id list for given may not have cds feature. If not I want to try again with a different id. May need to execute get_sequence from another function
# and pass id into it

	for feature in genbank_record.features:
			if feature.type == 'CDS':
				cds_feature = feature
				break
	if cds_feature:
		sequence = cds_feature.location.extract(genbank_record.seq)
		# print(description,sequence)
	else:

		print(f'cds_feature error: {organism_name}')
		return None

	return (description,str(sequence))


fc = fc.Fasta_cross()
fc.get_names(read_file = 'enterobacterales_short.fa',fasta_type = 'alex')

bacteria_names = fc.bacteria_names
gene_name = 'BamA'
records = [get_sequence(organism_name = bacteria_name,gene_name = gene_name) for bacteria_name in bacteria_names]
records = [record for record in records if record != None]
print(records)

with open('results.fasta','w') as f:
	for record in records:
		f.write(f'>{record[0]}\n{record[1]}\n')


# Only getting less than half back...