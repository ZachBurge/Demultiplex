#!/usr/bin/env python

import argparse
import bioinfo
import numpy as np
import matplotlib.pyplot as plt
import gzip

def get_args():
    '''Function for setting command line flags'''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="name of file used for input")
    parser.add_argument("-o", "--output", help="name of file used for output")
    parser.add_argument("-l", "--length", help="length of input file sequences")
    return parser.parse_args()

args = get_args()

# R1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
# R2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
# R3 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
# R4 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'

nucleotide_dists = np.zeros((int(args.length), 42), dtype=int)
print(nucleotide_dists)

with gzip.open(args.input, "rt") as infile:
    for i, line in enumerate(infile):
        if i%4 == 3:
            line = line.strip()
            for index, base in enumerate(line):
                q_score = bioinfo.convert_phred(base)
                nucleotide_dists[index, q_score] += 1

# print(nucleotide_dists)

phred_scores = np.arange(42) #[0,1,2,3,4,5,...,41]
total_scores_per_pos = np.sum(nucleotide_dists*phred_scores, axis=1)
total_reads = 363246735
mean = total_scores_per_pos/total_reads
# print(mean)

# mean = np.mean(nucleotide_dists, axis=1)
plt.scatter(range(int(args.length)), mean)
plt.title("Average Quality Score at each Nucleotide Position in Sequence Across all Reads")
plt.xlabel("Nucleotide")
plt.ylabel("Average Quality Score")
plt.savefig(f"{args.output}_hist.png")