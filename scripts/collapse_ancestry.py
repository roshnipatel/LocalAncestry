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
    parser.add_argument('--ind', help='Optional comma-delimited list of individual sample IDs for which to generate ancestry tract files. If not provided, ancestry tracts will be generated for all individuals in pop_map file.', default=None)
    parser.add_argument('--pop_labels', default='CEU,YRI', help='Comma-delimited list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels.')
    parser.add_argument('--out', help='Prefix to bed file. .{sample_ID}.A.bed and .{sample_ID}.B.bed will be appended.', required=True)
    parser.add_argument('--chr', help='Chromosome to process (should be an integer).', required=True)

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

def find_haplotype_bounds(indiv, index, add, pop_order, tracts, npop, chr):
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
                tracts.write(indiv + '\t' + str(chr) + '\t' + str(add) + '\t' + 'UNK' +
                            '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] +
                            '\t' + post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            else:
                tracts.write(indiv + '\t' + str(chr) + '\t' + str(add) + '\t' +
                            pop_order[int(last_anc_pos_cm[0])-1] + '\t' +
                            post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] + '\t' +
                            post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            post_anc_pos_cm = current_anc_pos_cm

        last_anc_pos_cm = current_anc_pos_cm

    # Last iteration, still need to print
    if last_anc_pos_cm[0] == -9:
        tracts.write(indiv + '\t' + str(chr) + '\t' + str(add) + '\t' + 'UNK' +
                    '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
                    '\t' + post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    else:
        tracts.write(indiv + '\t' + str(chr) + '\t' + str(add) + '\t' +
                    pop_order[int(current_anc_pos_cm[0])-1] + '\t' +
                    post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] + '\t' +
                    post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')

def main(current_ind, index, pop_order, chr):
    # Open bed files (2 haplotypes per individual)
    tracts = open(args.out, 'w')
    find_haplotype_bounds(current_ind, index, 0, pop_order, tracts, npop, chr)
    find_haplotype_bounds(current_ind, index, 1, pop_order, tracts, npop, chr)
    tracts.close()

if __name__ == '__main__':
    # Load parameters and files
    args = parse_args()
    pop_labels = args.pop_labels.split(',')
    npop = len(pop_labels)
    chr = int(args.chr)

    # Determine ordering of individuals in RFMix input
    ind_order = open(args.pop_map).readline().strip().split('\t')

    if args.ind:
        ind_list = args.ind.strip().split(',')
    else:
        ind_list = ind_order

    # Generate ancestry tract files for each individual in pop_map file
    for indiv in ind_list:
        idx = ind_order.index(indiv)
        print 'Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
        main(indiv, idx, pop_labels, chr)
        print 'Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
