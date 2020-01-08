"""
Estimates global ancestry proportion from local ancestry tracts.

Modified lightly from Alicia Martin (armartin via GitHub).
"""

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--bed', help='List of bed files.', nargs='+', required=True)
parser.add_argument('--ind', help='List of individual sample IDs for which to generate global ancestry estimates.', nargs='+', required=True)
parser.add_argument('--pops', default='EUR,AFR', help='Comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels', required=True)
parser.add_argument('--out', help='Filepath for output', required=True)

args = parser.parse_args()
out = open(args.out, 'w')
pops = args.pops.strip().split(',')
out.write('ID\t' + '\t'.join(pops) + '\n')
lai_props = [0]*len(pops)

# Extract individual sample ID from bed filepath
def extract_ID(filename):
    s = filename.split(".")
    ID = s[2]
    return(ID)

props = {}
for ID in args.ind:
    props[ID] = [0]*len(pops)

for bed_file in args.bed:
    ID = extract_ID(bed_file)
    if ID in args.ind:
        for tract in open(bed_file):
            tract = tract.strip().split()
            if tract[3] in pops: # This excludes stuff not listed in pops
                props[ID][pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))

for ID, anc_count in props.items():
    anc_prop = [round(i/sum(anc_count), 4) for i in anc_count]
    out.write(ID + '\t' + '\t'.join([str(i) for i in anc_prop]) + '\n')

out.close()
