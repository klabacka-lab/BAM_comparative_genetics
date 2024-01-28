import fasta_cross as fc
import re
from Bio import Entrez,SeqIO

Entrez.email = 'david.bean@utahtech.edu'


'''
Using fasta_cross to pull bacteria names from a fasta file
'''
def pull_names(read_file,fasta_type):
	fasta_cross = fc.Fasta_cross()
	fasta_cross.get_names(read_file = read_file,fasta_type = 'alex')
	return fasta_cross.bacteria_names


'''

'''
def build_query(gene_name,names_list,exclude_list):

	names_list = [f'"{name}"[Primary Organism]' for name in names_list]
	search_query = "("+' OR '.join(names_list)+f" AND {gene_name}[Gene])"
	if len(exclude_list) != 0:
		exclude_list = [f'"{name}"' for name in exclude_list]
		search_query = search_query + " NOT ("+" OR ".join(exclude_list)+")"

	# print('SEARCH QUERY: ',search_query)
	return search_query





def get_records(search_query,retmax):

	esearch_handle = Entrez.esearch(
		db ='protein',
		term = search_query,
		usehistory="y",
		retmax = retmax
		)

	e_record = Entrez.read(esearch_handle)
	esearch_handle.close()
	

	idlst = e_record['IdList']
	print(f'{len(idlst)} ids found')
	post = Entrez.epost('protein',id = idlst) 

	try:
		search_results = Entrez.read(post)
	except:
		'No results Error'
		return None
	webenv = search_results["WebEnv"]
	query_key = search_results["QueryKey"]


	fetch_handle =  Entrez.efetch(
		db = 'protein',
		rettype = "genbank",
		retmode = "fasta",
		webenv = webenv,
		query_key = query_key
		)

	records = list(SeqIO.parse(fetch_handle,"genbank"))

	return records




def sift_records(records,search_list,found_list,intruder_list):

	records_string = ""
	print('Intruder(s):',intruder_list)
	

	for record in records:
		for name in search_list:

			found_name = re.search(r'\[(.*)\]',record.description)
			try: 
				found_name = found_name.group(0)
				found_name = found_name[1:-1]# Snipping off the brackets []
				

			except:
				print("AttributeError: 'Nonetype' object has no attribute 'group;")
				break
			 

			if found_name in search_list and found_name not in found_list and found_name not in intruder_list:
				found_list.append(found_name)
				records_string += f'>{record.description}\n{record.seq}\n'
				break
			else:
				if found_name not in intruder_list:
					intruder_list.append(found_name)
				print("Intruder Found:", found_name)
				break

	print ("Intruder list: ",intruder_list)
				


	return (records_string,found_list,intruder_list)



def main():

	database = 'protein'
	gene_name = 'BamA'
	read_file = 'enterobacterales.fa'
	fasta_type = 'alex'

	names_list = pull_names(read_file,fasta_type)

	names_status = {name : False for name in names_list}
	search_list = [name for name in names_status if names_status[name] == False]

	results_string = ""
	searches = 0

	found_list = []
	intruder_list = []

	
	'''
 	Checking number of times the same ammount of results has been returned could be bad. It's possible to get 3 results back 3 times in a row even if the results aren't the same each time.
  	Need to check the actual content of the results.
   	I guess this is a band-aid for now
 	'''
	prev_number_names_found = 0
	same_number_found = 0
	

	write_file = open("search_results.fasta","w")

	while searches < len(names_list) and len(search_list) > 0 and same_number_found < 3:

		total_names = len(names_list)
		number_names_found = len(found_list)
		

		if prev_number_names_found == number_names_found:
			same_number_found += 1
		else:
			same_number_found = 0


		prev_number_names_found = number_names_found

		search_query = build_query(gene_name,search_list,exclude_list = intruder_list)
		records = get_records(search_query,retmax = 15)
		sift_results = sift_records(records,search_list,found_list,intruder_list)
		searches += 1

		# Writing found sequence to file
		write_file.write(sift_results[0])

		# removing found names from the search list
		for item in sift_results[1]:
			if item not in found_list:
				found_list.append(item)

		for name in sift_results[2]:
			if name not in intruder_list:
				intruder_list.append(name)

		print('Search Count: ',searches)
		for name in found_list:
			names_status[name] = True

		# Rebuilding search list without bacteria names already found
		search_list = [name for name in names_status if names_status[name] == False]

		print(f'{number_names_found} out of {total_names} names found\n')

	write_file.close()

	# Reporting names not found
	print ("\nFinished Searching\n")
	if len(search_list) != 0:
		missing_names = "\n".join(search_list)
		print (f'failed to find: {missing_names}\nCheck for outdate names. Writing missing names to missing_names.txt')
		with open('missing_names.txt','w') as file:
			file.write(missing_names)

	# Writing search results to a file
	
if __name__ == '__main__':
	main()





# https://biopython-tutorial.readthedocs.io/en/latest/notebooks/09%20-%20Accessing%20NCBIs%20Entrez%20databases.html#Using-the-history-and-WebEnv
# https://bioinformatics.stackexchange.com/questions/21788/getting-all-protein-entries-from-entrez-for-a-particular-taxonomic-id


# Maybe I could do only a few searches. If I get a certain species I could exclude it from the species 
# names list on the next pass. If I call 100 results at a time this way I may only have to do a few searches.


# Gotta keep a list of names that I get from search results that are not present in my list and explicitely exclude them from future
# searches
