"""
Combines estimated global ancestry fractions from RFMix and ADMIXTURE across all
individuals.
"""

import argparse
import pandas as pd

def calculate_error(row):
    avg = (row.YRI_admixture + row.YRI_rfmix) / 2
    error = ((row.YRI_admixture - avg) ** 2 + (row.YRI_rfmix - avg) ** 2) ** .5
    return(error)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix', help='File with global ancestry inferred from rfmix output. Should have columns ID, CEU, YRI.', required=True)
    parser.add_argument('--admixture', help='File with admixture output. Headerless, need to infer columns.', required=True)
    parser.add_argument('--map', help='Headerless, 2-row file with 1 column per individual: sample ID, race/ancestry classification. Individuals should be in same order as RFMix input.', required=True)
    parser.add_argument('--pops', help="Comma-delimited list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels.", required=True)
    parser.add_argument('--out', help="Filepath for output", required=True)
    args = parser.parse_args()
    
    # Split desired ancestries into a list
    pops = args.pops.strip().split(',')
    
    # Read file with RFMix global ancestry inference
    rfmix = {}
    with open(args.rfmix, 'r') as f:
        header = True
        r = f.readline().strip().split('\t')
        rfmix[r[0]] = r[1:]
        for row in f:
            r = row.strip().split('\t')
            rfmix[r[0]] = [float(x) for x in r[1:]]
    rfmix = pd.DataFrame(rfmix).set_index('ID').T
    
    # Read population map file
    with open(args.map, 'r') as f:
        sample_IDs = f.readline().strip().split('\t')
        sample_anc = f.readline().strip().split('\t')
    
    samples = list(zip(sample_IDs, sample_anc))
    
    # Create dictionary to store ADMIXTURE ancestry fractions for the reference
    # populations. This is because ADMIXTURE outputs ancestry fractions in a random
    # order, so we need to do some sleuthing (by looking at the reference populations)
    # to determine which column corresponds to which ancestry in the output.
    ref_pop_order = {k: [0] * len(pops) for k in pops}
    
    # Read file with ADMIXTURE global fractions
    admix = {}
    with open(args.admixture, 'r') as f:
        for ID, anc in samples:
            row = f.readline().strip().split()
            frac = [float(x) for x in row]
            admix[ID] = frac
            if anc in pops: # Checks for whether individual belongs to a reference population
                curr_idx = frac.index(max(frac)) # Determines largest ancestry proportion
                ref_pop_order[anc][curr_idx] += 1 # Documents largest ancestry proportion for this reference individual
    
    # For each reference population that RFMix/ADMIXTURE was run with, determine the
    # index (in ADMIXTURE output data) that best describes the individuals in this
    # population.
    order = []
    for pop in pops:
        max_idx = ref_pop_order[pop].index(max(ref_pop_order[pop]))
        order.append(max_idx)
    
    # Here we check to make sure that each reference population was assigned a different
    # column/index in the ADMIXTURE output. Otherwise, who knows what happened.
    if len(order) != len(set(order)):
        raise Exception('Reference populations did not have unambiguous assignment.')
    else:
        # Determine the correct ordering of population labels for the ADMIXTURE output
        admix['ID'] = list(zip(*sorted(list(zip(pops, order)), key=lambda x:x[1])))[0]
        admix = pd.DataFrame(admix).set_index('ID').T
    
    # Merge data by individual sample IDs
    comb = pd.merge(admix, rfmix, left_index=True, right_index=True, suffixes=('_admixture', '_rfmix'))
    
    # Identify individuals whose ADMIXTURE and RFMix calls are discordant
    comb["error"] = comb.apply(calculate_error, axis=1)
    
    comb.to_csv(args.out, sep='\t')
