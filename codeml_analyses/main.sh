#!/bin/bash

#SBATCH --time=150:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10240M   # memory per CPU core
#SBATCH -J "CODEML"   # job name
#SBATCH --mail-user=smithga@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


module load miniforge3
mamba activate paml_env

codeml 'control.ctl' 
