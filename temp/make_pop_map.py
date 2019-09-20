import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--samples', help='list of samples on which local and global ancestry inference will be performed. should be in same order as rfmix input.', nargs='+')
parser.add_argument('--admix', help='path to admixed sample data. requirements: file is csv, first column is indiv. sample ID.', required=True)
parser.add_argument('--admix_pop_col', help='name of column that specifies ancestry in admixed sample data', required=True)
parser.add_argument('--admix_pop_val', help='comma-delimited values of ancestry column in admixed data to filter for', required=True)
parser.add_argument('--ref', help='path to reference sample data csv. requirements: file is csv, first column is indiv. sample ID.', required=True)
parser.add_argument('--ref_pop_col', help='name of column that specifies ancestry in reference sample data', required=True)
parser.add_argument('--ref_pop_val', help='comma-delimited values of ancestry column in reference data to filter for', required=True)
parser.add_argument('--out', help='output file', required=True)
args = parser.parse_args()

sample_anc_map = {}

admix_pop_val = args.admix_pop_val.split(',')
ref_pop_val = args.ref_pop_val.split(',')

admix = pd.read_csv(args.admix, index_col=0)
admix_filtered = admix[admix[args.admix_pop_col].isin(admix_pop_val)]
for index, row in admix_filtered.iterrows():
    sample_anc_map[index] = 'ADMIX'

ref = pd.read_csv(args.ref, index_col=0)
ref_filtered = ref[ref[args.ref_pop_col].isin(ref_pop_val)]
for index, row in ref_filtered.iterrows():
    sample_anc_map[index] = row[args.ref_pop_col]

ancestries = []
for ID in args.samples:
    ancestries.append(sample_anc_map[ID])

anc_ID_map = {k: v + 1 for v, k in enumerate(ref_pop_val)}
anc_ID_map['ADMIX'] = 0

with open(args.out, 'w') as f:
    f.write("\t".join(args.samples))
    f.write('\n')
    f.write("\t".join(ancestries))
    f.write('\n')
    f.write("\t".join([anc_ID_map[k] for k in ancestries]))
