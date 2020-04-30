import sys
import pandas as pd
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', help='vcf with original statistical phasing', required=True)
parser.add_argument('--rfmix_phased', help='filepath for RFMix re-phased alleles', required=True)
parser.add_argument('--out', help='output filepath for differentially phased alleles', required=True)
args = parser.parse_args()

def check_phasing(row):
    if row.OG_phasing == row.RFMix_phasing:
        return(False)
    return(True)

geno_list = []
with gzip.open(args.vcf, 'rt') as geno_file:
    for line in geno_file:
        if line[0] == "#":
            continue
        else:
            line = line.split()
            pos = int(line[1])
            genotypes = line[9:]
            genotypes = "".join(genotypes).translate({ord('|'): ''})
            geno_list.append([pos, genotypes])
OG_phasing = pd.DataFrame(geno_list, columns=["Position", "OG_phasing"])

rfmix_phasing = pd.read_csv(args.rfmix_phased, names=["RFMix_phasing"])
phasing_df = pd.concat([OG_phasing, rfmix_phasing], axis=1)
phasing_df["Match"] = phasing_df.apply(check_phasing, axis=1)
phasing_df[phasing_df.Match == False].to_csv(args.out, sep='\t', index=False)
