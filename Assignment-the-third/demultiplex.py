#!/usr/bin/env python

import gzip
import itertools
import numpy as np
import bioinfo
import math
import argparse

def reverse_complement(seq: str):
    '''Function that returns the reverse complement of a given sequence'''
    complement = str.maketrans('ACGTN', 'TGCAN')
    return seq[::-1].translate(complement)

def meets_qscore_cutoff(seq: str, cutoff: float):
    '''Function that takes a sequence, in this case of phred quality scores,
        and returns True if each score in the sequence meets the given cutoff, False if not'''
    return all(bioinfo.convert_phred(b) >= cutoff for b in seq)

def write_record(file, record: list):
    '''Function that takes a file header and a list of values and writes the list to the file.'''
    file.writelines(f"{record}\n")

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cutoff", help="the quality score cutoff")
    return parser.parse_args()

args = get_args()

indexes_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
known_indexes = []
with open(indexes_file, "r") as fh:
    for line in fh:
        line = line.strip().split()
        known_indexes.append(line[4]) # the actual index sequence
known_indexes = known_indexes[1:] # remove the header line from the list

known_index_pairs = {}
output_files = {}
for i in itertools.product(known_indexes, repeat=2): # all possible pair combinations of the indexes
    pair = f"{i[0]}-{i[1]}"
    known_index_pairs[pair] = 0 # initialize each pair, matching and non matching
    if i[0] == i[1]: # when they are matching, need to create an output file
        output_files[f"R1-{pair}"] = open(f"outputs/R1-{pair}.fastq", "w")
        output_files[f"R2-{pair}"] = open(f"outputs/R2-{pair}.fastq", "w")
# create the rest of the output files
output_files["R1-hopped"] = open("outputs/R1-hopped.fastq", "w")
output_files["R2-hopped"] = open("outputs/R2-hopped.fastq", "w")
output_files["R1-unknown"] = open("outputs/R1-unknown.fastq", "w")
output_files["R2-unknown"] = open("outputs/R2-unknown.fastq", "w")

unknown = 0
# R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
# R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
R2 = "../TEST-input_FASTQ/test_R2.fastq"
R3 = "../TEST-input_FASTQ/test_R3.fastq"
R4 = "../TEST-input_FASTQ/test_R4.fastq"
R1 = "../TEST-input_FASTQ/test_R1.fastq"

# with gzip.open(R1, "rt") as R1, gzip.open(R2, "rt") as R2, gzip.open(R3, "rt") as R3, gzip.open(R4, "rt") as R4:
with open(R1, "rt") as R1, open(R2, "rt") as R2, open(R3, "rt") as R3, open(R4, "rt") as R4:
    while True:
        R1_header = R1.readline()
        if R1_header == '':
            break
        R1_seq = R1.readline()
        R1.readline() # the '+'
        R1_qscores = R1.readline()
        R1_lines = [R1_header, R1_seq, '+\n', R1_qscores]

        R2.readline() # header is same as R1, not needed
        R2_seq = R2.readline()
        R2.readline() # the '+'
        R2_qscores = R2.readline()
        R2_lines = [R1_header, R2_seq, '+\n', R2_qscores]

        R3.readline() # header is same as R4, not needed
        R3_seq = R3.readline()
        R3.readline() # the '+'
        R3_qscores = R3.readline()

        R4_header = R4.readline()
        R4_seq = R4.readline()
        R4.readline() # the '+'
        R4_qscores = R4.readline()
        R3_lines = [R4_header, R3_seq, '+\n', R3_qscores]
        R4_lines = [R4_header, R4_seq, '+\n', R4_qscores]

        R3_seq = R3_seq.strip()
        print(R3_seq)
        R3_rev_comp = reverse_complement(R3_seq)
        print(R3_rev_comp)

        if f"{R2_seq}-{R3_rev_comp}" not in known_index_pairs:
            output_files["R1-unknown"].writelines(R1_lines)
            output_files["R2-unknown"].writelines(R4_lines)

for k, v in output_files.items():
    v.close()