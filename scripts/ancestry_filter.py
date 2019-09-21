"""
Reads sample metadata to identify sample IDs corresponding to the specified ancestry.
"""

import sys, csv
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sample_data', help='Path to sample metadata. Requirements: file is csv, first column is indiv. sample ID in VCF.', required=True)
parser.add_argument('--pop_col', help='Name of column that specifies race/ancestry in admixed sample data.', required=True)
parser.add_argument('--pop_val', help='Comma-delimited values of race/ancestry column to filter for.', required=True)
parser.add_argument('--out', help='Filepath of output file.', required=True)
args = parser.parse_args()

# Split desired ancestries into a list
pop_val = args.pop_val.strip().split(',')

# Reads in sample metadata, setting first column as index. Indices of filtered
# rows will be written to output file, which is why sample IDs must be the first
# column in the metadata file.
samples = pd.read_csv(args.sample_data, index_col=0)
filtered_samples = samples[samples[args.pop_col].isin(pop_val)]
filtered_IDs = set(list(filtered_samples.index.values))

# Write sample IDs of individuals with specified race/ancestry to file.
with open(args.out, 'w') as out:
    for ID in filtered_IDs:
        out.write(ID)
        out.write("\n")
