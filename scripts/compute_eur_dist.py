"""
Parses through bed files for all chromosomes to compute empirical distribution of European tract lengths.
"""

import argparse
import pandas as pd
import os
from lai_global import extract_ID

parser = argparse.ArgumentParser()
parser.add_argument('--bed', help='Master directory containing bed files for each chromosome.', required=True)
parser.add_argument('--tract_out', help="Filepath for tract length output", required=True)
args = parser.parse_args()

all_bed_files = []
for i in range(1, 23):
    chr_path = os.path.join(args.bed, "chr{0}".format(str(i)))
    chr_files = [os.path.join(chr_path, f) for f in os.listdir(chr_path) if f[-3:] == "bed"]
    all_bed_files.extend(chr_files)

eur_tract_lengths = pd.DataFrame()
for bf in all_bed_files:
    bf_df = pd.read_csv(bf, sep='\t', names=['chrom', 'gpos_start', 'gpos_stop', 'anc', 'pos_start', 'pos_stop'])
    idv = extract_ID(bf)
    eur_tracts = bf_df[bf_df.anc == "CEU"]
    if not eur_tracts.empty:
        tmp = eur_tracts.apply(lambda row: [idv, row.pos_stop - row.pos_start], axis=1, result_type='expand')
        eur_tract_lengths = pd.concat([eur_tract_lengths, tmp])
eur_tract_lengths.columns = ["nwd_id", "length"]

eur_tract_lengths.to_csv(args.tract_out, sep='\t', index=False)

