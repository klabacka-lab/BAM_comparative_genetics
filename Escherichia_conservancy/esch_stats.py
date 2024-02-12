
import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
import pandas as pd
import os, sys
import warnings





# Creating result data frame
new_df = pd.DataFrame({})



# Iterating through files in stats folder

os.chdir('./stats')
files = os.listdir()

totalFiles = len(files)
count = 0
print( f'{count}/{totalFiles} genes processed ({round((count/totalFiles)*100)}%)')


for fileName in files:

	geneData = {}

	# Obtaining Gene Name
	gene_name = fileName.split("_")[0]
	geneData['GENE'] = gene_name

	df = pd.read_csv(fileName)

	# Obtaining mean Proportion
	# Rounding to 4th decimal arbitrarily. Can change later
	meanProportion = df[["Proportion"]].mean().iloc[0].round(4)
	medProportion = df[["Proportion"]].median().iloc[0]
	geneData['Mean Prop'] = [meanProportion]
	geneData['Med Prop'] = [medProportion]


	# Merging Gene info to outfile Data Frame
	geneData = pd.DataFrame(geneData)


	if new_df.empty:
		new_df = geneData
	else:
		new_df = pd.concat([new_df,geneData],ignore_index=True)

	count+=1

	sys.stdout.write("\033[F")
	print(f'{count}/{totalFiles} genes processed ({(count/totalFiles)*100}%)')

# Returning to original working directory and writing results csv
os.chdir('..')
new_df.to_csv('Gene_Proportions.csv',sep=',',index=False)
