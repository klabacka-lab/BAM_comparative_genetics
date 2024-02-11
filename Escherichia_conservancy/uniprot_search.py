import requests
import random
import fasta_cross
import sys
import time

# Pulls gene names from ecoli names file. Returns names as list
def get_genes(inFile):
    with open(inFile,'r') as file:
        geneList = [i.split()[0] for i in file.readlines()]
        return(geneList)

# requesting FASTA sequences from uniprot
def uniprotSearch(gene,taxID):
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A{taxID}%29+AND+%28gene%3A{gene}%29%29'
    return requests.get(url).text

def write_searches(gene):
    for gene in randomSample:
        with open(f'./search_results/{gene}_enterobacterales.fasta','w') as outFile:
            outFile.write(uniprotSearch(gene))

# I think this is limiting the number of Escherichia species I can retrieve because ultimately it will
# only return species that were already in enterobacterales.fa. This late in the game we are going to just
# run with the list of bacteria we have though.

def cross_fasta(names_fasta,search_fasta,write_handle):
    crosser = fasta_cross.Fasta_cross()
    crosser.get_names(names_fasta,'alex') # enterobacterales.fa
    crosser.cross_ref(search_fasta,crosser.bacteria_names,repeat=False)
    crosser.write_fasta(crosser.match_records,handle = write_handle)
    # We will algn sequences in a later step with align.py

# def cross_fasta(search_fasta,bactList):
#     crosser = fasta_cross.Fasta_cross()
#     crosser.cross_ref(search_fasta,bactList,handle=write_handle)

def main():

    taxID = sys.argv[2]
    geneHandle = sys.argv[1]

    # Limiting number of genes retrieved for debugging
    if sys.argv[3]:
        maxGenes = int(sys.argv[3]) + 1
    else:
        maxGenes = None

    if maxGenes == None:
        geneList = get_genes(geneHandle)
    else:
        geneList = get_genes(geneHandle)[1:maxGenes:]

    print(f'Retrieving genes from {geneHandle}')
    print(f'Max genes: {maxGenes - 1}')


    totalGenes = len(geneList)
    progressCount = 0
    failureCount = 0
    startTime = time.time()

    print(f'{progressCount}/{totalGenes} searched')

    for gene in geneList:
        try:
            write_handle = f'./search_results/{gene}_enterobacterales.fasta'
            searchResult = uniprotSearch(gene,taxID) 
            with open(write_handle,'w') as outFile:
                    outFile.write(searchResult)

            cross_fasta(
                    names_fasta = 'enterobacterales.fa',
                    search_fasta = write_handle, 
                    write_handle = write_handle
                    )

            progressCount += 1
            currentTime = time.time()
            elapsedTime = round(currentTime - startTime)

            sys.stdout.write("\033[F")
            print(f'{progressCount}/{totalGenes} searched. Time elapsed: {elapsedTime} s. Search Failures {failureCount}')
        except:
            failureCount += 1
            print(f'Search Failures: {failureCount}')

if __name__ == '__main__':
    main()

