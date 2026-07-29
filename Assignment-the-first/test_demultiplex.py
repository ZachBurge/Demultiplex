#!/usr/bin/env python

import itertools

indexes = ['TCTTCGAC','AACAGCGA','GTCCTAAG','TATGGCAC','TACCGGAT','ACGATCAG','TCGACAAG','ATCCGGTA']
indexes_reversed = indexes[::-1]

known_index_pairs = {}
for j in itertools.product(indexes, repeat=2):
    known_index_pairs[f"{j[0]}-{j[1]}"] = 0
    
