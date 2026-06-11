This directory contains the images generated from pymol for the proteins we analysed.

The color_residues.py python scripts contain the code to generate images, and can be ran after the pymol_env conda environment is activated.

To create the pymol_env conda environment, use this command in terminal:
```
conda create -n pymol_env python=3.11 pymol-open-source -c conda-forge
```

To activate the environment, use this command in terminal:
```
conda activate pymol_env
```

You can then run the script using this code:
```
python color_residues.py
```
