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
    return seq.translate(complement[::-1])

def meets_qscore_cutoff(seq: str, cutoff: float):
    '''Function that takes a sequence, in this case of phred quality scores,
        and returns True if each score in the sequence meets the given cutoff, False if not'''
    return all(bioinfo.convert_phred(b) >= cutoff for b in seq)

def write_record(file, record: list):
    '''Function that takes a file header and a list of values and writes the list to the file.'''
    pass