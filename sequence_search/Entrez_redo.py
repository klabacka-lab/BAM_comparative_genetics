from Bio import SeqIO,Entrez
Entrez.email = "david.bean@utahtech.edu"
database = "gene"
#term = "txid562[Organism] AND b0321[Gene] AND biomol_mrna[PROP]"
term = "txid562[Organism] AND bamA[Gene]"
searchHandle = Entrez.esearch(
	db = database,
	retmax = 10,
	term = term
	)

searchRecord = Entrez.read(searchHandle)
searchHandle.close()
idList = searchRecord['IdList']



fetchHandle = Entrez.efetch(
	db = 'nucleotide',
	id = idList[0],
	rettype="fasta",
	retmode = "text"
	)


record = fetchHandle.read()
print(record)

