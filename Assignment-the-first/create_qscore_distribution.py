#!/usr/bin/env python

import argparse
import bioinfo
import numpy as np

def get_args():
    '''Function for setting command line flags'''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="name of file used for input")
    parser.add_argument("-o", "--output", help="name of file used for output")
    return parser.parse_args()

args = get_args()

# R1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
# R2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
# R3 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
# R4 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'

nucleotide_dists = np.zeros((101, 363246735), dtype=float)

# with open(args.input, "r") as infile, open(args.output, "w") as outfile:
