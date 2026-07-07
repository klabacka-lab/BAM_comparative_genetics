#!/usr/bin/env python3
"""
Map FEL codon positions (alignment columns) to 8ADG CIF residue numbers.

The FEL codon positions refer to columns in the amino acid alignment
(codon_align.faa). This script counts ungapped E. coli residues up to each
column, then adds the 3-residue MGS-tag offset to get the CIF residue number.

Usage:
  python map_codons.py codon_positions.txt
  python map_codons.py codon_positions.txt --alignment ../HyPhy/codon_align.faa
  python map_codons.py codon_positions.txt -o cif_residues.txt
"""

import argparse
from pathlib import Path

from Bio import SeqIO


ALIGNMENT  = Path(__file__).parent.parent / 'HyPhy' / 'codon_align.faa'
ECOLI_ID   = 'Escherichia_coli'
CIF_OFFSET = 3   # MGS tag at N-terminus of 8ADG reference
CIF_FIRST  = 27  # first residue resolved in the CIF structure
CIF_LAST   = 813 # last residue resolved in the CIF structure


def build_mapping(alignment_path: Path) -> dict:
    """Return {alignment_col (1-based): cif_residue or None}."""
    ecoli_seq = None
    for rec in SeqIO.parse(alignment_path, 'fasta'):
        if rec.id == ECOLI_ID:
            ecoli_seq = str(rec.seq)
            break

    if ecoli_seq is None:
        raise ValueError(f"'{ECOLI_ID}' not found in {alignment_path}")

    mapping = {}
    ungapped = 0

    for col, aa in enumerate(ecoli_seq, start=1):
        if aa != '-':
            ungapped += 1
            cif_res = ungapped + CIF_OFFSET
            mapping[col] = cif_res if CIF_FIRST <= cif_res <= CIF_LAST else None
        else:
            mapping[col] = None

    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', type=Path,
                        help='Codon positions file (comma-separated or one per line)')
    parser.add_argument('--alignment', type=Path, default=ALIGNMENT,
                        help=f'Amino acid alignment FASTA (default: {ALIGNMENT})')
    parser.add_argument('-o', '--output', type=Path, default=Path('output/cif_residues.txt'),
                        help='Output file for CIF residue numbers (default: output/cif_residues.txt)')
    args = parser.parse_args()

    mapping = build_mapping(args.alignment)

    codon_positions = []
    with open(args.input) as fh:
        for token in fh.read().replace('\n', ',').split(','):
            token = token.strip()
            if token:
                codon_positions.append(int(token))

    cif_residues = []
    skipped = []
    for pos in codon_positions:
        cif = mapping.get(pos)
        if cif is not None:
            cif_residues.append(cif)
        else:
            skipped.append(pos)

    with open(args.output, 'w') as fh:
        for r in cif_residues:
            fh.write(f"{r}\n")

    print(f"Mapped {len(cif_residues)} codon positions to CIF residues -> {args.output}")
    if skipped:
        print(f"Skipped {len(skipped)} positions (gap in E. coli or outside resolved structure): "
              f"{skipped[:10]}{'...' if len(skipped) > 10 else ''}")


if __name__ == '__main__':
    main()
