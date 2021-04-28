"""
Generates RFMix v1.5.4 input files for p and q arms of a single chromosome. See
RFMix v1.5.4 manual for details on file formats.
"""

import sys, os
import pandas as pd
from genetic_map import map_positions
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--admix', help='Path to admixed sample metadata. Requirements: file is csv, first column is indiv. sample ID.', required=True)
parser.add_argument('--admix_id_col', help='Name of column that specifies VCF sample ID in admixed sample metadata.', required=True)
parser.add_argument('--admix_pop_col', help='Name of column that specifies race/ancestry in admixed sample metadata.', required=True)
parser.add_argument('--admix_pop_val', help='Comma-delimited values of race/ancestry column to filter for.', required=True)
parser.add_argument('--ref', help='Path to reference sample metadata. Requirements: file is csv, first column is indiv. sample ID.', required=True)
parser.add_argument('--ref_id_col', help='Name of column that specifies VCF sample ID in reference sample metadata.', required=True)
parser.add_argument('--ref_pop_col', help='Name of column that specifies race/ancestry in reference sample metadata.', required=True)
parser.add_argument('--ref_pop_val', help='Comma-delimited values of race/ancestry column in reference data to filter for.', required=True)
parser.add_argument('--genetic_map', help='Path to map from genetic positions to chromosomal coordinates.', required=True)
parser.add_argument('--out', help='Prefix for output files.', required=True)
parser.add_argument('--cent', help='Starting position of centromere for chromosome being processed.', required=True)
args = parser.parse_args()

admix_pop_val = args.admix_pop_val.split(',')
ref_pop_val = args.ref_pop_val.split(',')

cent_pos = int(args.cent)

# Parse input VCF and write alleles file. Store ordered lists of SNP positions
# and sample IDs.
with open(args.out + ".alleles.p.txt", 'w') as p_file, open(args.out + ".alleles.q.txt", 'w') as q_file:
    p_pos = []
    q_pos = []
    for line in sys.stdin:
        if line[:2] == "##":
            continue
        if line[0] == "#":
            line = line.split()
            samples = line[9:]
        else:
            line = line.split()
            genotypes = line[9:]
            genotypes = "".join(genotypes).translate({ord('|'): ''})
            pos = int(line[1])
            if pos < cent_pos: # Determine which arm of chromosome SNP is located on
                p_file.write(genotypes + "\n")
                p_pos.append(pos)
            else:
                q_file.write(genotypes + "\n")
                q_pos.append(pos)

# Create snp locations file (LF-delimited file that contains genetic position of
# every marker). Also write a map from SNP chromosomal position to genetic position.
genetic_map = pd.read_csv(args.genetic_map, sep='\s+', names=["Chrom", "SNP", "GeneticDist", "Pos"])

q_gen_pos = map_positions(genetic_map, q_pos)
with open(args.out + ".snp_locations.q.txt", 'w') as snp_loc:
    for gp in q_gen_pos:
        snp_loc.write(str(gp))
        snp_loc.write('\n')
with open(args.out + ".pos_map.q.txt", "w") as f:
    f.write("\t".join([str(p) for p in q_pos]))
    f.write("\n")
    f.write("\t".join([str(p) for p in q_gen_pos]))

if p_pos: # Chromosome is not acrocentric, and both arms can be analyzed in parallel
    p_gen_pos = map_positions(genetic_map, p_pos)
    with open(args.out + ".snp_locations.p.txt", 'w') as snp_loc:
        for gp in p_gen_pos:
            snp_loc.write(str(gp))
            snp_loc.write('\n')
    with open(args.out + ".pos_map.p.txt", "w") as f:
        f.write("\t".join([str(p) for p in p_pos]))
        f.write("\n")
        f.write("\t".join([str(p) for p in p_gen_pos]))
else: # Chromosome is acrocentric, and we only generate one set of output files
    os.remove(args.out + ".alleles.p.txt")

# Create classes file (space-delimited file that designates ancestry of every
# haplotype). Also write a map from sample IDs to ancestry.
sample_anc_map = {}

anc_ID_map = {k: v + 1 for v, k in enumerate(ref_pop_val)}
anc_ID_map['ADMIX'] = 0

admix = pd.read_csv(args.admix, index_col=args.admix_id_col, sep='\t')
admix_filtered = admix[admix[args.admix_pop_col].isin(admix_pop_val)]
for index, row in admix_filtered.iterrows():
    sample_anc_map[index] = 'ADMIX'

ref = pd.read_csv(args.ref, index_col=args.ref_id_col, sep='\t')
ref_filtered = ref[ref[args.ref_pop_col].isin(ref_pop_val)]
for index, row in ref_filtered.iterrows():
    sample_anc_map[index] = row[args.ref_pop_col]

ancestries = []
ancestry_IDs = []
for s in samples:
    ancestries.append(sample_anc_map[s])
    ancestry_IDs.append(anc_ID_map[sample_anc_map[s]])

with open(args.out + ".classes.txt", 'w') as f:
    classes = [val for val in ancestry_IDs for _ in (0, 1)]
    f.write(' '.join([str(a) for a in classes]))

with open(args.out + ".pop_map.txt", 'w') as f:
    f.write("\t".join(samples))
    f.write('\n')
    f.write("\t".join(ancestries))
    f.write('\n')
    f.write("\t".join([str(i) for i in ancestry_IDs]))
