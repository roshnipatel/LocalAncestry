"""
Estimates global ancestry proportion from local ancestry tracts.

Modified lightly from Alicia Martin (armartin via GitHub).
"""

import argparse
import os

def extract_ID(filename):
    ID = filename.split('/')[-1].split('.')[0]
    return(ID)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed', help='Directory containing bed files.', required=True)
    parser.add_argument('--pops', help='Comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels', required=True)
    parser.add_argument('--out', help='Filepath for output', required=True)
    args = parser.parse_args()
    
    out = open(args.out, 'w')
    pops = args.pops.strip().split(',')
    out.write('ID\t' + '\t'.join(pops) + '\n')
    
    
    all_bed_files = [os.path.join(args.bed, f) for f in os.listdir(args.bed) if f[-3:] == "bed"]
    
    props = {}
    for bed_file in all_bed_files:
        ID = extract_ID(bed_file) # Extract individual sample ID from bed filepath
        for tract in open(bed_file):
            tract = tract.strip().split()
            if tract[3] in pops: # This excludes stuff not listed in pops
                if ID in props:
                    props[ID][pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))
                else:
                    props[ID] = [0]*len(pops)
                    props[ID][pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))
    
    for ID, anc_count in props.items():
        anc_prop = [round(i/sum(anc_count), 4) for i in anc_count]
        out.write(ID + '\t' + '\t'.join([str(i) for i in anc_prop]) + '\n')
    
    out.close()
