#!/usr/bin/env python3
"""
Parse and filter HyPhy FEL output to extract significant codon positions.

Columns in the FEL table:
  Codon | Partition | alpha | beta | LRT | Selection detected?

At least one primary filter (--pvalue or --omega) is required.

Usage examples:
  # Filter by p-value
  python parse_fel.py ../HyPhy/HyPhy_FEL_output.txt --pvalue 0.05

  # Filter by dN/dS (omega = beta/alpha); keeps sites with omega < value
  python parse_fel.py ../HyPhy/HyPhy_FEL_output.txt --omega 0.5

  # Combine both
  python parse_fel.py ../HyPhy/HyPhy_FEL_output.txt --pvalue 0.05 --omega 0.5

  # Only negative selection sites
  python parse_fel.py ../HyPhy/HyPhy_FEL_output.txt --omega 0.5 --selection Neg

  # Launch interactive Streamlit visualization
  python parse_fel.py ../HyPhy/HyPhy_FEL_output.txt --vis
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


_HEADER_RE = re.compile(r'\|\s*Codon\s*\|')
_ROW_RE = re.compile(
    r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
    r'\s*(Neg|Pos)\.\s*p\s*=\s*([\d.]+)\s*\|'
)


def parse_fel_output(filepath: Path) -> pd.DataFrame:
    """Parse HyPhy FEL output file and return a DataFrame of significant sites."""
    records = []
    in_table = False
    with open(filepath) as fh:
        for line in fh:
            if _HEADER_RE.search(line):
                in_table = True
                continue
            if not in_table:
                continue
            m = _ROW_RE.search(line)
            if m:
                records.append({
                    'Codon':     int(m.group(1)),
                    'Partition': int(m.group(2)),
                    'alpha':     float(m.group(3)),
                    'beta':      float(m.group(4)),
                    'LRT':       float(m.group(5)),
                    'Selection': m.group(6),
                    'pvalue':    float(m.group(7)),
                })
    return pd.DataFrame(records)


def apply_filters(df: pd.DataFrame, args) -> pd.DataFrame:
    mask = pd.Series(True, index=df.index)
    if args.pvalue is not None:
        mask &= df['pvalue'] <= args.pvalue
    if args.omega is not None:
        with np.errstate(invalid='ignore', divide='ignore'):
            omega = np.where(df['alpha'] > 0, df['beta'] / df['alpha'], np.nan)
        mask &= np.isfinite(omega) & (omega < args.omega)
    if args.partition is not None:
        mask &= df['Partition'] == args.partition
    if args.alpha_min is not None:
        mask &= df['alpha'] >= args.alpha_min
    if args.alpha_max is not None:
        mask &= df['alpha'] <= args.alpha_max
    if args.beta_min is not None:
        mask &= df['beta'] >= args.beta_min
    if args.beta_max is not None:
        mask &= df['beta'] <= args.beta_max
    if args.lrt_min is not None:
        mask &= df['LRT'] >= args.lrt_min
    if args.lrt_max is not None:
        mask &= df['LRT'] <= args.lrt_max
    if args.selection != 'both':
        mask &= df['Selection'] == args.selection
    return df[mask]


def write_codon_positions(codons: list, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as fh:
        fh.write(','.join(str(c) for c in sorted(codons)) + '\n')
    print(f"Wrote {len(codons)} codon positions to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Parse and filter HyPhy FEL output to extract codon positions.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('input', type=Path,
                        help='Path to HyPhy FEL output file')
    parser.add_argument('-o', '--output', type=Path, default=Path('output/codon_positions.txt'),
                        help='Output file for filtered codon positions (default: output/codon_positions.txt)')

    filter_group = parser.add_argument_group('filters (at least --pvalue or --omega is required)')
    filter_group.add_argument('--pvalue', type=float, metavar='FLOAT',
                              help='Maximum p-value (e.g. 0.05)')
    filter_group.add_argument('--omega', type=float, metavar='FLOAT',
                              help='Maximum dN/dS (ω = β/α); keeps sites with ω < value. '
                                   'Sites with α = 0 (undefined ω) are excluded.')
    filter_group.add_argument('--selection', choices=['Neg', 'Pos', 'both'], default='both',
                              help='Selection type: Neg, Pos, or both (default: both)')
    filter_group.add_argument('--partition', type=int, metavar='INT',
                              help='Keep only sites in this partition number')
    filter_group.add_argument('--alpha-min', type=float, metavar='FLOAT',
                              help='Minimum alpha (synonymous substitution rate)')
    filter_group.add_argument('--alpha-max', type=float, metavar='FLOAT',
                              help='Maximum alpha')
    filter_group.add_argument('--beta-min', type=float, metavar='FLOAT',
                              help='Minimum beta (non-synonymous substitution rate)')
    filter_group.add_argument('--beta-max', type=float, metavar='FLOAT',
                              help='Maximum beta')
    filter_group.add_argument('--lrt-min', type=float, metavar='FLOAT',
                              help='Minimum LRT (likelihood ratio test statistic)')
    filter_group.add_argument('--lrt-max', type=float, metavar='FLOAT',
                              help='Maximum LRT')

    parser.add_argument('--vis', action='store_true',
                        help='Launch interactive Streamlit visualization with cutoff sliders')

    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: file not found: {args.input}")

    if not args.vis and args.pvalue is None and args.omega is None:
        sys.exit("Error: at least one of --pvalue or --omega is required.\n"
                 "  Examples:\n"
                 "    make wrangle PVALUE=0.05\n"
                 "    make wrangle OMEGA=0.5\n"
                 "    make wrangle PVALUE=0.05 OMEGA=0.5")

    df = parse_fel_output(args.input)
    print(f"Parsed {len(df)} significant sites from {args.input.name}")

    if args.vis:
        vis_app = Path(__file__).parent / 'fel_vis_app.py'
        try:
            subprocess.run(
                [sys.executable, '-m', 'streamlit', 'run', str(vis_app),
                 '--', str(args.input.resolve())],
                check=True,
            )
        except FileNotFoundError:
            sys.exit("Error: streamlit not found. Install it with: pip install streamlit")
        return

    filtered = apply_filters(df, args)
    print(f"After filtering: {len(filtered)} sites retained")

    if filtered.empty:
        print("Warning: no sites matched the filter criteria.")

    write_codon_positions(filtered['Codon'].tolist(), args.output)


if __name__ == '__main__':
    main()
