import os,sys,re
from Bio import SeqIO

# {'Escherichia albertii': 1757, 'Escherichia coli': 2730, 'Escherichia fergusonii': 2037, 'Escherichia marmotae': 2102, 'Escherichia ruysiae': 2131}

def main():
    path = "./search_results"
    resultsDir = os.listdir(path)
    specRegex = 'OS=([A-Za-z]* [a-z]*)'

    speciesCount = {} # Dictionary {(key) species : (value) # occurances among fasta files)

    genesOfInterest = []
    speciesOfInterest = [
            'Escherichia albertii',
            'Escherichia coli',
            'Escherichia fergusonii',
            'Escherichia ruysiae'
            ]

    for fileHandle in resultsDir:
        filePath = path+"/"+fileHandle
        with open(filePath,'r') as f:

            # Counting occurrances of each species across all files
            fasta_species = re.findall(specRegex,f.read())
            for species in fasta_species:
                    if species in speciesCount.keys():
                        speciesCount[species] = speciesCount[species] + 1
                    else:
                        speciesCount[species] = 1

            # Recording genes that have all 5 species of interest
            match = True
            for species in speciesOfInterest:
                    if species not in fasta_species:
                        match = False
            if match == True:
                    genesOfInterest.append(fileHandle) # appending file Handle if it contais genes from all species of interest 
            sys.stdout.write("\033[F")
            print(speciesCount)
    print('copying genes of interest')
    
    filterCount = 0
    minLength = 100
    for handle in genesOfInterest:
        unfilteredPath = f'./search_results/{handle}'
        filteredPath = f'./filtered_results/{handle}'
        records = list(SeqIO.parse(unfilteredPath,'fasta'))

        tooShort = False

        for record in records:
            if len(record.seq) < minLength:
                tooShort = True
                print(handle + 'too short')

        if tooShort == False:
            filterCount += 1
            os.system(f'cp {unfilteredPath} {filteredPath}')
    print(str(filterCount) + ' genes remaining after filtering for min length ' + str(minLength))
if __name__ == '__main__':
    main()
