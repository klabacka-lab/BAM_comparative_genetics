import os
import fasta_cross as fc
fastas = os.listdir("./filtered_results")

crosser = fc.Fasta_cross()

for fasta in fastas:
	crosser.align(f'./filtered_results/{fasta}',f'./filtered_results/{fasta}')
	