import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--rfmix', help='file with global ancestry inferred from rfmix output. should have columns ID, CEU, YRI.', required=True)
parser.add_argument('--admixture', help='file with admixture output. headerless, need to infer columns.', required=True)
parser.add_argument('--map', help='pop map file.', required=True)
parser.add_argument('--pops', help="order of populations in rfmix output.", required=True)
parser.add_argument('--out', help="output file", required=True)
args = parser.parse_args()

pops = args.pops.strip().split(',')

rfmix = {}
with open(args.rfmix, 'r') as f:
    header = True
    for row in f:
        r = row.strip().split('\t')
        rfmix[r[0]] = r[1:]
rfmix = pd.DataFrame(rfmix).set_index('ID').T
print(rfmix)

with open(args.map, 'r') as f:
    sample_IDs = f.readline().strip().split('\t')
    sample_anc = f.readline().strip().split('\t')

samples = list(zip(sample_IDs, sample_anc))

ref_pop_order = {k: [0] * len(pops) for k in pops} # Accumulate evidence for whether ADMIXTURE outputted the correct population ordering
admix = {}
with open(args.admixture, 'r') as f:
    for ID, anc in samples:
        row = f.readline().strip().split()
        frac = [float(x) for x in row]
        admix[ID] = frac
        if anc in pops:
            curr_idx = frac.index(max(frac))
            ref_pop_order[anc][curr_idx] += 1
print(admix)

order = []
for pop in pops:
    max_idx = ref_pop_order[pop].index(max(ref_pop_order[pop]))
    order.append(max_idx)
print(order)

if len(order) != len(set(order)):
    raise Exception('Reference populations did not have unambiguous assignment.')
else:
    admix['ID'] = list(zip(*sorted(list(zip(pops, order)), key=lambda x:x[1])))[0]
    admix = pd.DataFrame(admix).set_index('ID').T

comb = pd.merge(admix, rfmix, left_index=True, right_index=True, suffixes=('_admixture', '_rfmix'))
comb.to_csv(args.out, sep='\t')
