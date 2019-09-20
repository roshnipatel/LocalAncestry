import sys, csv
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sample_data', help='path to sample data. requirements: file is csv, first column is indiv. sample ID.', required=True)
parser.add_argument('--pop_col', help='name of column that specifies ancestry in admixed sample data', required=True)
parser.add_argument('--pop_val', help='comma-delimited values of ancestry column to filter for', required=True)
parser.add_argument('--out', help='output file', required=True)
args = parser.parse_args()

# Parse sample metadata to extract IDs of individuals with desired ancestry
pop_val = args.pop_val.strip().split(',')
samples = pd.read_csv(args.sample_data, index_col=0)
filtered_samples = samples[samples[args.pop_col].isin(pop_val)]
filtered_IDs = set(list(filtered_samples.index.values))

with open(args.out, 'w') as out:
    for ID in filtered_IDs:
        out.write(ID)
        out.write("\n")
