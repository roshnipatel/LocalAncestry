"""
Combines estimated global ancestry fractions from RFMix and ADMIXTURE across all
chromosomes.
"""

from combine_ancestry_calls import calculate_error
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--anc', help='Files with global ancestry inferred from RFMix and ADMIXTURE.', required=True, nargs='+')
parser.add_argument('--map', help='Filepath for information on chromosome length.', required=True)
parser.add_argument('--frac_out', help="Filepath for global ancestry fraction output", required=True)
parser.add_argument('--idv_out', help="Filepath for individual exclusion list", required=True)
args = parser.parse_args()

filepath = args.anc[0]
anc = pd.read_csv(filepath, delimiter='\t', index_col=0).reset_index().rename({"index": "ID"}, axis=1)
chrom = int(filepath.strip().split('/')[1][3:])
anc['Chrom'] = chrom

for filepath in args.anc[1:]:
    tmp = pd.read_csv(filepath, delimiter='\t', index_col=0).reset_index().rename({"index": "ID"}, axis=1)
    chrom = int(filepath.strip().split('/')[1][3:])
    tmp['Chrom'] = chrom
    anc = pd.concat([anc, tmp])
anc = anc.drop("error", axis=1)

# Multiply per-chromosome ancestry fractions by chromosome lengths and normalize
# whole genome sum to 1 (to obtain whole genome ancestry fractions)
chr_map = pd.read_csv(args.map, delimiter='\t')
merged = pd.merge(anc, chr_map, left_on='Chrom', right_on='Chrom')
merged["YRI_admixture"] = merged["YRI_admixture"] * merged["ChromosomeEnd"]
merged["CEU_admixture"] = merged["CEU_admixture"] * merged["ChromosomeEnd"]
merged["YRI_rfmix"] = merged["YRI_rfmix"] * merged["ChromosomeEnd"]
merged["CEU_rfmix"] = merged["CEU_rfmix"] * merged["ChromosomeEnd"]
total_bp = chr_map['ChromosomeEnd'].sum()
global_frac = merged.groupby('ID').sum().divide(total_bp)

# Identify individuals to exclude from downstream analyses due to high 
# discordance between RFMix and ADMIXTURE global ancestry estimates
global_frac["error"] = global_frac.apply(calculate_error, axis=1)
exclusion_idv = set(global_frac[global_frac.error > 0.05].index)
exclusion_idv = exclusion_idv.union(set(global_frac[global_frac.YRI_admixture < 0.5].index))
with open(args.idv_out, 'a+') as f:
    for idv in exclusion_idv:
        f.write(idv)
        f.write("\n")

global_frac[["YRI_admixture", "CEU_admixture", "YRI_rfmix", "CEU_rfmix", "error"]].to_csv(args.frac_out, sep='\t', float_format='%.5f', )
