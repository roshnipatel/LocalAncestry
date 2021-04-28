"""
Takes in RFMix Viterbi output and generates ancestry tracts for each haplotype
of each individual specified in input.

Modified lightly from Alicia Martin (armartin via GitHub).
"""

from datetime import datetime
import time
import argparse
from itertools import izip_longest
import gzip
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix', help='Path to RFMix Viterbi output, chr expected in filename.', required=True)
    parser.add_argument('--snp_map', help='Headerless, 2-row file with 1 column per SNP: bp, cM. Should only contain SNPs that RFMix was run with, and SNPs should be in same order as RFMix input.', required=True)
    parser.add_argument('--fbk', default=None)
    parser.add_argument('--fbk_threshold', type=float, default = 0.9)
    parser.add_argument('--pop_map', help='Headerless, 2-row file with 1 column per individual: sample ID, race/ancestry classification. Individuals should be in same order as RFMix input.', required=True)
    parser.add_argument('--pop_labels', default='CEU,YRI', help='Comma-delimited list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels.')
    parser.add_argument('--out', help='Prefix to bed file. {sample_ID}.A.bed and {sample_ID}.B.bed will be appended.', required=True)
    parser.add_argument('--chrom', help='Chromosome to process (should be an integer).', required=True)

    args = parser.parse_args()
    return(args)

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def check_gt_posterior(fbk_max, fbk_threshold, index, add, hap_anc, line, current_anc):
    if fbk_max[index*2+add] >= fbk_threshold:
      hap_anc = rfmix_viterbi[index*2+add]
    else:
      hap_anc = -9
      current_anc[add] = -9
    return (current_anc, hap_anc)

def find_haplotype_bounds(index, add, pop_order, hap, npop, chr):
    print str(chr) + ' [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    rfmix_file = re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.rfmix)
    if rfmix_file.endswith('gz'):
        rfmix = gzip.open(rfmix_file)
    else:
        rfmix = open(rfmix_file)
    if args.fbk is not None:
        fbk_file = re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.fbk)
        if fbk_file.endswith('gz'):
            fbk = gzip.open(fbk_file)
        else:
            fbk = open(fbk_file)
    snp_map = open(re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.snp_map))
    bp = snp_map.readline().strip().split('\t')
    cM = snp_map.readline().strip().split('\t')

    last_anc_pos_cm = [None, None, 0]

    counter = 0
    for line in rfmix:
        cur_pos, cur_gpos = bp[counter], cM[counter]
        counter += 1
        rfmix_viterbi = line.strip().split()
        if args.fbk is not None:
            fbk_line = fbk.readline().strip().split()
            fbk_line = map(float, fbk_line)
            fbk_max = []
            for i in grouper(npop, fbk_line):
                fbk_max.append(max(i))
            if fbk_max[index*2+add] < args.fbk_threshold:
                rfmix_viterbi[index*2+add] = -9
        # Fencepost for start of the chromosome
        if counter == 1:
            last_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
            post_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
            continue

        # Start regular iterations
        current_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
        if current_anc_pos_cm[0] == last_anc_pos_cm[0]:
            last_anc_pos_cm = current_anc_pos_cm
            continue
        else:
            # We've reached the end of a region. Need to print.
            if last_anc_pos_cm[0] == -9:
                hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] +
                        '\tUNK\t' + post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            else:
                hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] + '\t' +
                        pop_order[int(last_anc_pos_cm[0])-1] + '\t' +
                        post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            post_anc_pos_cm = current_anc_pos_cm

        last_anc_pos_cm = current_anc_pos_cm

    # Last iteration, still need to print
    if last_anc_pos_cm[0] == -9:
        hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
            '\tUNK\t' + post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    else:
        hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
            '\t' + pop_order[int(current_anc_pos_cm[0])-1] + '\t' +
            post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')

def main(current_ind, index, pop_order, chr):
    # Open bed files (2 haplotypes per individual)
    hap_a = open(args.out + current_ind + '.A.bed', 'w')
    hap_b = open(args.out + current_ind + '.B.bed', 'w')

    find_haplotype_bounds(index, 0, pop_order, hap_a, npop, chr)
    hap_a.close()

    find_haplotype_bounds(index, 1, pop_order, hap_b, npop, chr)
    hap_b.close()

if __name__ == '__main__':
    # Load parameters and files
    args = parse_args()
    pop_labels = args.pop_labels.split(',')
    npop = len(pop_labels)
    chrom = int(args.chrom)

    # Identify admixed individual IDs for which to generate bed files
    with open(args.pop_map, 'r') as f:
        ind_ids = f.readline().strip().split('\t')
        pop_type = f.readline().strip().split('\t')

    admix_idv = {}
    i = 0
    for tup in zip(ind_ids, pop_type):
        if tup[1] != "ADMIX":
            break
        admix_idv[i] = tup[0]
        i += 1

    # Generate ancestry tract files for each individual specified in command-line argument
    for idx, ind_id in admix_idv.items():
        print 'Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
        main(ind_id, idx, pop_labels, chrom)
        print 'Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    status_file = open(args.out + "job_complete.txt", 'w')
    status_file.close()
