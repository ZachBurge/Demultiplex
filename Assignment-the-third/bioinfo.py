#!/usr/bin/env python

# Author: Zach Burge

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

import pytest

DNA_bases = {"A", "T", "G", "C", "N"}
RNA_bases = {"A", "U", "G", "C", "N"}

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """This function calculates the average Phred quality score of the given phred_score.
    It will iterate through each encoded score in the phred_score string and add its corresponding decoded score to the sum,
    before dividing that sum by the length of the phred_score string to get the average of the decoded phred scores."""
    sum = 0
    for score in phred_score:
        sum += convert_phred(score)
    return sum / len(phred_score)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    seq_set = set(seq)
    return RNA_bases >= seq_set if RNAflag else DNA_bases >= seq_set

def gc_content(seq, RNAflag=False):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    if validate_base_seq(seq, RNAflag=False):
        return (seq.count("G") + seq.count("C"))/len(seq)
    raise AssertionError

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    if len(lst)%2 != 0:
        return lst[len(lst)//2]
    else:
        return (lst[int((len(lst)/2)-1)]+lst[int((len(lst)/2))])/2

def oneline_fasta(file, newfile):
    '''docstring'''
    with open(file, "r") as infile, open(newfile, "w") as outfile:
        dna_line = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if dna_line != '':
                    outfile.write(f"{dna_line}\n")
                outfile.write(f"{line}\n")
                dna_line = ''
            else:
                dna_line += line
        outfile.write(dna_line)

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert qual_score('JJJJJJJJ') == 41, "qual_score is not returning the correct average"
    assert qual_score('GGGGGAAA') == 35.75, "qual_score is not returning the correct average"
    assert qual_score('CCFFFFFHHHHGIJIHHHIJJJJJIIHJHGGIBGGHIJGHHEIIJIHGIBIIIIIIJJJIGGHIJJIIIHGHGIGHIAHDE?ACE<@BBDDCD@5:<<C>') == 37.38, "qual_score is not returning the correct average"
    assert qual_score('CCFFFFDHFHHHJGJIGE?HIEHJBG)?DHEEFHIEIIIJJJJGGEIFHECCCDE;') == 36.732142857142854, "qual_score is not returning the correct average"
    print("Your qual_score function is working correctly")
    assert validate_base_seq('JJJJJJJJ') == False, "validate_base_seq is not working correctly"
    assert validate_base_seq('GGGGGAAA') == True, "validate_base_seq is not working correctly"
    assert validate_base_seq('CCFFFFFHHHHGIJIHHHIJJJJJIIHJHGGIBGGHIJGHHEIIJIHGIBIIIIIIJJJIGGHIJJIIIHGHGIGHIAHDE?ACE<@BBDDCD@5:<<C>') == False, "validate_base_seq is not working correctly"
    assert validate_base_seq('TESTINGTESTINGTESTINGTESTING') == False, "validate_base_seq is not working correctly"
    print("Your validate_base_seq function is working correctly")
    with pytest.raises(AssertionError):
        gc_content('JJJJJJJJ')
    assert gc_content('GGGGGAAA') == 0.625, "gc_content is not working correctly"
    with pytest.raises(AssertionError):
        gc_content('CCFFFFFHHHHGIJIHHHIJJJJJIIHJHGGIBGGHIJGHHEIIJIHGIBIIIIIIJJJIGGHIJJIIIHGHGIGHIAHDE?ACE<@BBDDCD@5:<<C>')
    assert gc_content('TCTTCTAATTTATCAAGCAATACTTATAAATCTTATCATCACACTCCTTTTTTAAGAGTA') == 0.25, "gc_content is not working correctly"
    print("Your gc_content function is working correctly")
    assert calc_median([4,4,3,4,5,6,7]) == 4, "calc_median does not work for odd length list"
    assert calc_median([1,2,50]) == 2, "calc_median does not work for odd length list"
    assert calc_median([4,3,4,5,6,7]) == 4.5, "calc_median does not work for even length list"
    assert calc_median([1,2,3,4,5,6,7,8,9,10]) == 5.5, "calc_median does not work for odd length list"
    print("Your calc_median function is working correctly")
