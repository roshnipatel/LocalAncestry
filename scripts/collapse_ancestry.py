__author__ = 'armartin'
#takes in RFMix output and creates collapsed 2 haploid bed files per individual

from datetime import datetime
import time
import argparse
from itertools import izip_longest
import gzip
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix', help='path to RFMix Viterbi output, chr expected in filename', required=True)
    parser.add_argument('--snp_map', help='headerless, 2-row file with 1 column per SNP: bp, cM', required=True)
    parser.add_argument('--fbk', default=None)
    parser.add_argument('--fbk_threshold', type=float, default = 0.9)
    parser.add_argument('--ind', help='space-delimited list of individual sample IDs', required=True, nargs='+')
    parser.add_argument('--pop_map', help='path to file mapping Individual IDs to ancestries', required=True)
    parser.add_argument('--pop_labels', default='ASN,EUR,AFR',
                    help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
    parser.add_argument('--out', help='prefix to bed file, .{sample_ID}.A.bed and .{sample_ID}.B.bed will be appended', required=True)
    parser.add_argument('--chr', help='chromosome to process (should be an integer)', default='', required=False)

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
        #fencepost for start of the chromosome
        if counter == 1:
            last_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
            post_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
            continue

        #start regular iterations
        current_anc_pos_cm = [rfmix_viterbi[2*index + add], cur_gpos, cur_pos]
        if current_anc_pos_cm[0] == last_anc_pos_cm[0]:
            last_anc_pos_cm = current_anc_pos_cm
            continue
        else:
            #we've reached the end of a region. Need to print.
            if last_anc_pos_cm[0] == -9:
                hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] +
                        '\tUNK\t' + post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            else:
                hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] + '\t' +
                        pop_order[int(last_anc_pos_cm[0])-1] + '\t' +
                        post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            post_anc_pos_cm = current_anc_pos_cm

        last_anc_pos_cm = current_anc_pos_cm

    #last iteration, still need to print
    if last_anc_pos_cm[0] == -9:
        hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
            '\tUNK\t' + post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    else:
        hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
            '\t' + pop_order[int(current_anc_pos_cm[0])-1] + '\t' +
            post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')

def main(current_ind, index, pop_order, chr):

    #open bed files (2 haplotypes per individual)
    hap_a = open(args.out + '.' + current_ind + '.A.bed', 'w')
    hap_b = open(args.out + '.' + current_ind + '.B.bed', 'w')

    find_haplotype_bounds(index, 0, pop_order, hap_a, npop, chr)
    hap_a.close()

    find_haplotype_bounds(index, 1, pop_order, hap_b, npop, chr)
    hap_b.close()

if __name__ == '__main__':
    #load parameters and files
    args = parse_args()
    pop_labels = args.pop_labels.split(',')
    npop = len(pop_labels)
    chr = int(args.chr)
    ind_order = open(args.pop_map).readline().strip().split('\t')
    for i in args.ind:
        idx = ind_order.index(i) # Grab index of current individual
        print 'Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
        main(i, idx, pop_labels, chr)
        print 'Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
