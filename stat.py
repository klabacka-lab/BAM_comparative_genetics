from position import Position

aa_dict = {}

with open("enterobacteralesOA.fa", "r") as fasta_file:
            read_file = fasta_file.read()
            lines = read_file.split(">")
            for line in lines:
                temp = line.split("\n")
                print(temp[0])
                print(temp[1])


