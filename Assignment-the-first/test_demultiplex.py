#!/usr/bin/env python

import itertools

indexes = ['TCTTCGAC','AACAGCGA','GTCCTAAG','TATGGCAC','TACCGGAT','ACGATCAG','TCGACAAG','ATCCGGTA']

for j in itertools.combinations_with_replacement(indexes, 2):
    print(f"{j[0]}-{j[1]}")

# known_indexes = {'TCTTCGAC-TCTTCGAC': 0, 'AACAGCGA-AACAGCGA': 0, 'GTCCTAAG-GTCCTAAG': 0, 'TATGGCAC-TATGGCAC': 0}