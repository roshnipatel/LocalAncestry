import sys, csv
import pandas as pd
from genetic_map import map_positions
import numpy as np

cent_file = sys.argv[1]
cent_table = pd.read_csv(cent_file, delimiter='\t', skiprows=1, names=["region", "chr", "start", "stop", "len"])
cent_table = cent_table[cent_table["region"].str.startswith('CEN')] # Ignore the region called HET7

chrom_file = sys.argv[2]
chrom_table = pd.read_csv(chrom_file, names=["chr", "length", "accession1", "accession2"])

# Use Plink-generated genetic maps to map centromere coordinates from base pairs
# to centiMorgans. Calculate the center of the centromere as the average of the
# two positions.
coord = []
for i in range(1, 23):
    genetic_map = "data/maps/plink.chr{0}.GRCh38.map".format(str(i))
    map = pd.read_csv(genetic_map, delimiter=' ', names=["Chrom", "SNP", "GeneticDist", "Pos"])
    cent_pos = list(cent_table[cent_table.chr == str(i)][["start", "stop"]].iloc[0])
    chrom_length = int(chrom_table[chrom_table.chr == str(i)]["length"].iloc[0].replace(",", ""))
    positions = [0] + cent_pos + [chrom_length]
    genetic_pos = map_positions(map, positions)
    start, end = max(0, genetic_pos[0]), genetic_pos[3]
    center = max(0, np.mean(genetic_pos[1:3]))
    if abs(center - start) < 2:
        coord.append([i, start, end])
    coord.append([i, start, center, end])

with open("data/maps/centromere.map", 'w') as f:
    writer = csv.writer(f)
    for chrom in coord:
        writer.writerow(chrom)
