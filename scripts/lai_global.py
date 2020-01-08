"""
Estimates global ancestry proportion from local ancestry tracts.

Modified lightly from Alicia Martin (armartin via GitHub).
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tracts', help='File of called local ancestry tracts for all individuals', required=True)
parser.add_argument('--pops', default='EUR,AFR', help='Comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels', required=True)
parser.add_argument('--out', help='Filepath for output', required=True)
args = parser.parse_args()

out = open(args.out, 'w')
pops = args.pops.strip().split(',')
out.write('ID\t' + '\t'.join(pops) + '\n')

props = {}
with open(args.tracts, 'r') as tract_file:
    line = tract_file.readline().strip.split('\t')
    if tract[3] in pops: # This excludes stuff not listed in pops
        ID = tract[0]
        if ID not in props:
            props[ID] = [0] * len(pops)
        props[ID][pops.index(tract[3])] += (float(tract[7]) - float(tract[6]))

for ID, anc_count in props.items():
    anc_prop = [round(i/sum(anc_count), 4) for i in anc_count]
    out.write(ID + '\t' + anc_prop[0] + '\t' + anc_prop[1] + '\n')

out.close()
