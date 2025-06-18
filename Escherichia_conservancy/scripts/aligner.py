import os,sys
import fasta_cross as fc
# NOTE: must have muscle installed


def main():
	# Obtaining names of fasta files in filtered results directory
	fastas = os.listdir("./filtered_results")

	# Initiating crosser object, which contains alignment method
	crosser = fc.Fasta_cross()

	# Tracking progress
	totalFiles = len(fastas)
	count = 0

	print(f'{count}/{totalFiles} alignments complete')
	for fasta in fastas:
		crosser.align(f'./filtered_results/{fasta}',f'./filtered_results/{fasta}',quiet = True)
		count +=1
		sys.stdout.write("\033[F")
		print(f'{count}/{totalFiles} alignments complete ({round((count/totalFiles)*100)}%)')

if __name__ == '__main__':
	main()