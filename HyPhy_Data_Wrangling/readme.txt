
==Brain Storm==
Goals
  - Parse ../HyPhy/HyPhy_FEL_output.txt
  - Wrangle to filter for codons of interest based on features.
    (Unsure which feature. Pretty sure it's P value in that last column)
  - Return Codon Positions of selected codons. (First column)

Bonus
  - Interactive visualization with sliders for cutoff values of feature set
    I know that is Mckinley's job, but what I have in mind might look good on
    my personal website.

==Plan==
  - Input:
      Path to HyPhy_FEL_output.txt

      Flags for value cutoffs for each of the 5 non-target featrues (Partition, alpha, beta, LRT,
      Selection detected. If flag is omitted, no filtering will be performed based on that feature

      flag to launch vis with imbeded streamlit-imbedded pymol session. Cutoff sliders for features
      option to save figure as a png with codon positions in a seperate file. This option should
      not require any flags for cutoff values.

  - Wrangling: Python Pandas simple filter statements
  - Primary Output: A text file with a list of codon positions from HyPhy_FEL_output.txt


