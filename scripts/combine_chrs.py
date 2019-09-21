"""
Combines estimated global ancestry fractions from RFMix and ADMIXTURE across all
chromosomes.
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--anc', help='Files with global ancestry inferred from RFMix and ADMIXTURE.', required=True, nargs='+')
parser.add_argument('--map', help='Filepath for information on chromosome length.', required=True)
parser.add_argument('--out', help="Filepath for output", required=True)
args = parser.parse_args()

filepath = args.anc[0]
anc = pd.read_csv(filepath, delimiter='\t')
anc.columns = ['ID', 'YRI_admixture', 'CEU_admixture', 'CEU_rfmix', 'YRI_rfmix']
chrom = filepath.strip().split('.')[1][3:]
anc['Chrom'] = chrom

for filepath in args.anc[1]:
    tmp = pd.read_csv(filepath, delimiter='\t')
    tmp.columns = ['ID', 'YRI_admixture', 'CEU_admixture', 'CEU_rfmix', 'YRI_rfmix']
    chrom = filepath.strip().split('.')[1][3:]
    tmp['Chrom'] = chrom
    anc = pd.concat([anc, tmp])

# Multiply per-chromosome ancestry fractions by chromosome lengths and normalize
# whole genome sum to 1 (to obtain whole genome ancestry fractions)
chr_map = pd.read_csv(args.map, delimiter='\t')
merged = pd.merge(anc, chr_map, left_on='Chrom', right_on='Chrom')
merged["YRI_admix_bp"] = merged["YRI_admixture"] * merged["ChromosomeEnd"]
merged["CEU_admix_bp"] = merged["CEU_admixture"] * merged["ChromosomeEnd"]
merged["YRI_rfmix_bp"] = merged["YRI_rfmix"] * merged["ChromosomeEnd"]
merged["CEU_rfmix_bp"] = merged["CEU_rfmix"] * merged["ChromosomeEnd"]
total_bp = chr_map['ChromsomeEnd'].sum()
global_frac = merged.groupby('ID').sum().divide(total_bp)

global_frac.to_csv(args.out, delimiter='\t', index_col=False)
