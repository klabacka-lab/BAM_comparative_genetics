import re, os
from sys import argv
from Bio import SeqIO

class Fasta_cross():
	def __init__(self):
		self.bacteria_names = []
		self.match_records = {}

	def get_names(self,read_file,fasta_type):

		regex = {
		'txt' : None,
		'uniprot' : r'\[(.*)\]',
		'genbank' : r'OS=(\S* \S*)'
		}
		
		if regex[fasta_type] == None:
			with open(read_file) as f:
				self.bacteria_names = [line.strip() for line in f.readlines()]

		with open(read_file) as f:
			bacteria_names = re.findall(regex[fasta_type],f.read())
			self.bacteria_names = bacteria_names
	
	def cross_ref(self,search_fasta,names,repeat = False):
		match_species = []
		records = SeqIO.parse(search_fasta,"fasta")
		records = {record.description:str(record.seq) for record in records}

		if repeat == True:
			for name in self.bacteria_names:
				for description in records.keys():
					if re.search(name,description):
						match_species.append(name)
						self.match_records[description] = records[description]

		if repeat == False:
			for name in self.bacteria_names:
				for description in records.keys():
					if name not in match_species and re.search(name,description):
						match_species.append(name)
						self.match_records[description] = records[description]

	def write_fasta(self,records_dict,handle = 'DefaultName.fasta'):
		with open(handle,'w') as f:
			for description,sequence in zip(records_dict.keys(),records_dict.values()):
				f.write('>'+description+'\n')
				f.write(sequence+'\n')

	def align(self,read_fasta):
		os.system(f"muscle -in {read_fasta} -out aligned_{read_fasta}")


if __name__ == '__main__':

	for i,argument in enumerate(argv):
		if argument == '-sort':
			sort_fasta = argv[i+1]
		if argument == '-names':
			names_fasta = argv[i+1]
		if argument == '-out':
			out_handle = argv[i+1]
		if argument == '-type':
			fasta_type = argv[i+1]

	obj = Fasta_cross()
	obj.get_names(names_fasta,fasta_type)
	obj.cross_ref(sort_fasta,obj.bacteria_names,repeat = False)
	obj.write_fasta(obj.match_records, handle = out_handle)
	obj.align(read_fasta = out_handle)


