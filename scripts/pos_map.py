import sys, os
import pandas as pd
from genetic_map import map_positions
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--genetic_map', help='path to map from genetic positions to chromosomal coordinates', required=True)
parser.add_argument('--out', help='prefix for output files', required=True)
args = parser.parse_args()

pos = []
for line in sys.stdin:
    if line[:2] == "##":
        continue
    if line[0] == "#":
        continue
    else:
        line = line.split()
        p = int(line[1])
        pos.append(p)

# Create snp locations file (LF-delimited file that contains genetic position of
# every marker). Also write a map from SNP chromosomal position to genetic position.
genetic_map = pd.read_csv(args.genetic_map, delimiter=' ', names=["Chrom", "SNP", "GeneticDist", "Pos"])

gen_pos = map_positions(genetic_map, pos)
with open(args.out + ".0.pos_map.merged.txt", "w") as f:
    f.write("\t".join([str(p) for p in pos]))
    f.write("\n")
    f.write("\t".join([str(p) for p in gen_pos]))
