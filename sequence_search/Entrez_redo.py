from Bio import SeqIO,Entrez
Entrez.email = "david.bean@utahtech.edu"

term = "BamA[Gene Name]"

searchHandle = Entrez.esearch(
	db = "nucleotide",
	retmax = 10,
	term = term
	)

searchRecord = Entrez.read(searchHandle)
searchHandle.close()
idList = searchRecord['IdList']



fetchHandle = Entrez.efetch(
	db = "nucleotide",
	id = idList[0],
	rettype="genbank",
	retmode = "text"
	)

record = SeqIO.read(fetchHandle,'genbank')
print(record.seq)

