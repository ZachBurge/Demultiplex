#!/usr/bin/env python

import gzip
import itertools
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np

RC_TABLE = str.maketrans('ACGTN', 'TGCAN') # use global translation table so that it doesn't recreate it every time the function is called
def reverse_complement(seq: str):
    '''Function that returns the reverse complement of a given sequence'''
    return seq[::-1].translate(RC_TABLE)

def meets_qscore_cutoff(seq: str, cutoff: float):
    '''Function that takes a sequence, in this case of phred quality scores,
        and returns True if each score in the sequence meets the given cutoff, False if not'''
    # choosing not to use bioinfo.convert_phred to save on having to call that function in every iteration of the loop
    limit = cutoff + 33
    return all(ord(b) >= limit for b in seq)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cutoff", help="the quality score cutoff")
    return parser.parse_args()

args = get_args()
cutoff = int(args.cutoff)

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
        output_files[f"R1-{pair}"] = open(f"/scratch/bgmp/zburge/demux/R1-{pair}.fastq", "w")
        output_files[f"R4-{pair}"] = open(f"/scratch/bgmp/zburge/demux/R4-{pair}.fastq", "w")
# create the hopped and unknown output files
output_files["R1-hopped"] = open("/scratch/bgmp/zburge/demux/R1-hopped.fastq", "w")
output_files["R4-hopped"] = open("/scratch/bgmp/zburge/demux/R4-hopped.fastq", "w")
output_files["R1-unknown"] = open("/scratch/bgmp/zburge/demux/R1-unknown.fastq", "w")
output_files["R4-unknown"] = open("/scratch/bgmp/zburge/demux/R4-unknown.fastq", "w")

# caching the output file dictionary lookups to reference within the loop, saves a little time
R1_hopped = output_files["R1-hopped"]
R4_hopped = output_files["R4-hopped"]
R1_unknown = output_files["R1-unknown"]
R4_unknown = output_files["R4-unknown"]

unknown = 0
R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
# R1 = "../TEST-input_FASTQ/test_R1.fastq"
# R2 = "../TEST-input_FASTQ/test_R2.fastq"
# R3 = "../TEST-input_FASTQ/test_R3.fastq"
# R4 = "../TEST-input_FASTQ/test_R4.fastq"

record_num = 0 # counter to track how many records have been read
with gzip.open(R1, "rt") as R1, gzip.open(R2, "rt") as R2, gzip.open(R3, "rt") as R3, gzip.open(R4, "rt") as R4:
# with open(R1, "rt") as R1, open(R2, "rt") as R2, open(R3, "rt") as R3, open(R4, "rt") as R4:
    while True:
        R1_header = R1.readline().strip()
        if R1_header == '':
            break
        R1_seq = R1.readline()
        R1.readline() # the '+'
        R1_qscores = R1.readline()

        R2.readline() # header is same as R1, not needed
        R2_seq = R2.readline().strip()
        R2.readline() # the '+'
        R2_qscores = R2.readline().strip()

        R3.readline() # header is same as R4, not needed
        R3_seq = R3.readline().strip()
        R3.readline() # the '+'
        R3_qscores = R3.readline().strip()

        R4_header = R4.readline().strip()
        R4_seq = R4.readline()
        R4.readline() # the '+'
        R4_qscores = R4.readline()

        # to track progress in slurm.out
        record_num += 1
        if record_num % 100000 == 0:
            print(f"Processed {record_num} records", flush=True)

        R3_rev_comp = reverse_complement(R3_seq)

        R1_header += f" {R2_seq}-{R3_rev_comp}\n" # add index labels to header
        R4_header += f" {R2_seq}-{R3_rev_comp}\n" # add index labels to header
        R1_lines = [R1_header, R1_seq, '+\n', R1_qscores] # create list of lines to write later
        R4_lines = [R4_header, R4_seq, '+\n', R4_qscores] # create list of lines to write later

        pair = R2_seq + "-" + R3_rev_comp # used for dictionary lookup throughout loop
        if pair not in known_index_pairs: # index pair is unknown
            unknown += 1
            R1_unknown.writelines(R1_lines)
            R4_unknown.writelines(R4_lines)
        else: # the index pair is known, but might not be high enough quality
            if not meets_qscore_cutoff(R2_qscores, cutoff) or not meets_qscore_cutoff(R3_qscores, cutoff):
                # too low quality = write to unknown
                unknown += 1
                R1_unknown.writelines(R1_lines)
                R4_unknown.writelines(R4_lines)
            else:
                known_index_pairs[pair] += 1
                if R2_seq == R3_rev_comp: # indexes match, write to appropriate output file
                    output_files[f"R1-{pair}"].writelines(R1_lines)
                    output_files[f"R4-{pair}"].writelines(R4_lines)
                else: # indexes both exist, but do not match, write to hopped output file
                    R1_hopped.writelines(R1_lines)
                    R4_hopped.writelines(R4_lines)

# close all output files
for k, v in output_files.items():
    v.close()

# report all useful counts
matched = hopped = 0
heatmap = np.zeros((len(known_indexes), len(known_indexes)), dtype=int) # create empty 2D array for heatmap axes
indexes_to_heatmap = {index: i for i, index in enumerate(known_indexes)} # initialize empty dict for getting index coordinates on heatmap
with open("demux_report.txt", "w") as out:
    out.write("Index Pair\tCount\n")
    for pair, count in known_index_pairs.items():
        out.write(f"{pair}\t{count}\n")
        r2, r3 = pair.split('-')
        row_val = indexes_to_heatmap[r2] # x coordinate of heatmap
        col_val = indexes_to_heatmap[r3] # y coordinate of heatmap
        heatmap[row_val, col_val] = count # heatmap at x,y is the count of that index pairing
        if r2 == r3:
            matched += count
        else:
            hopped += count
    out.write(f"Total Matched: {matched}\n")
    out.write(F"Total Hopped: {hopped}\n")
    out.write(f"Total Unknown: {unknown}\n")

plt.figure(figsize=(10,8)) # slightly bigger than default
sns.heatmap(heatmap, norm=LogNorm(), xticklabels=known_indexes, yticklabels=known_indexes, cmap="viridis") # LogNorm because the values on the diagonal are much higher than any others, so without it the graph colors would look weird
plt.xlabel("Forward indexes")
plt.ylabel("Reverse indexes (rev comp)")
plt.title("Matched and Hopped Index Pair Counts")

plt.savefig("index_counts_heatmap.png", dpi=300) # dpi = pixels per inch, default is 100