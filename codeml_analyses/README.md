# **CodeML Analyses**

This part of the repository was created using the CODEML function of the PAML program.
They are organized as such:

# Nomenclature
Each file is named according to its lineage model (M) and site model (o). The numbers are associated with internal assignments of the models. 
For example, m0o0 is the file for constant estimated omega across both sites and lineages. 

# Contents
Each file contains several things:

- Slurm batch info: these files, starting with slurm-(integer) are the console readouts of the file. This info is mostly useful for debugging.
- Control files: these files, starting with the model followed by control.ctl, are the settings for each of the readouts.
- Out.txt: The fully compiled outputs of the program.
- other files: other files comprise parts of the out.txt, used in subanalyses.

# Other resources
For a further explanation, see the following manual:
https://academic.oup.com/mbe/article/40/4/msad041/7140562?login=true

If, for some reason, that link breaks, here is the DOI:
https://doi.org/10.1093/molbev/msad041
