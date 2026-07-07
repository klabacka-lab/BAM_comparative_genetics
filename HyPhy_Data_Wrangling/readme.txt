
==Overview==
HyPhy FEL data wrangling pipeline for identifying evolutionarily conserved
residues in BamA. Filters HyPhy FEL output by p-value and/or dN/dS (omega),
maps selected codon positions to 8ADG CIF residue numbers, and provides an
interactive Streamlit visualization with a 3D protein viewer.

The key insight: omega (dN/dS = beta/alpha) is likely the most biologically
meaningful filter. Sites with low omega are under strong purifying selection
-- the amino acid is intolerant of change across 99 bacterial lineages and
is therefore critical for BamA structure or function.


==Project Status==
Minimum requirements met. Pending integration testing.


==Files==
  parse_fel.py        -- Parse and filter HyPhy_FEL_output.txt
  map_codons.py       -- Map filtered codon positions to 8ADG CIF residue numbers
  fel_vis_app.py      -- Streamlit visualization app
  hist_component/     -- Custom Streamlit component: interactive dN/dS histogram
  Makefile            -- Targets: install, wrangle, map, vis, clean
  requirements.txt    -- Python dependencies (managed via .venv)


==Quick Start==
  make wrangle OMEGA=0.5          # filter by dN/dS < 0.5
  make wrangle PVALUE=0.05        # filter by p-value < 0.05
  make wrangle PVALUE=0.05 OMEGA=0.5  # combine both filters
  make map                        # map codon positions to CIF residue numbers
  make vis                        # launch Streamlit app


==Inputs==
  ../HyPhy/HyPhy_FEL_output.txt  -- HyPhy FEL results (pre-computed)
  ../HyPhy/codon_align.faa        -- Codon alignment used for the FEL run
  ../StructuralModeling/8adg.cif  -- Crystal structure for 3D visualization


==Outputs==
  output/codon_positions.txt  -- Comma-separated codon positions passing filters
  output/cif_residues.txt     -- CIF residue numbers for structure visualization


==Background: What is HyPhy FEL?==
HyPhy FEL (Fixed Effects Likelihood) was run externally on a codon alignment
of 99 bacterial BamA sequences (1,191 codon positions). For each site, HyPhy
estimated two substitution rates across the full phylogenetic tree:

  alpha -- synonymous rate (silent DNA changes, e.g. AAA->AAG, both Lys)
  beta  -- non-synonymous rate (amino acid-changing, e.g. AAA->GAA, Lys->Glu)

These per-site alpha and beta values are stored in HyPhy_FEL_output.txt.
The HyPhy output only includes sites where beta != alpha was statistically
significant (LRT p <= 0.1). This pipeline reads those values and computes:

  omega = beta / alpha  (dN/dS)

Sites with alpha = 0 have undefined omega and are excluded from omega filtering.


==Makefile Variables==
  INPUT       Path to FEL output (default: ../HyPhy/HyPhy_FEL_output.txt)
  PVALUE      Maximum p-value cutoff (e.g. PVALUE=0.05)
  OMEGA       Maximum dN/dS cutoff (e.g. OMEGA=0.5); keeps sites with omega < value
  SELECTION   Neg, Pos, or both (default: both)
  ALPHA_MIN / ALPHA_MAX   Synonymous rate bounds
  BETA_MIN  / BETA_MAX    Non-synonymous rate bounds
  LRT_MIN   / LRT_MAX     Likelihood ratio test statistic bounds
  PARTITION               Restrict to a specific partition number
  CODON_FILE  Output path for codon positions (default: output/codon_positions.txt)
  CIF_FILE    Output path for CIF residues (default: output/cif_residues.txt)
