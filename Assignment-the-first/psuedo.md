Functions needed:
```
def reverse_complement(seq: str): -> rev_comp_str: str
    '''Function that returns the reverse complement of a given sequence'''
    For each char in seq, convert it to its complement. If the char is an N, leave it as an N. After each char is complemented, reverse the string with slicing.
    return rev_seq
Input: 'AATAA'
Expected Output: 'TTATT'

def meets_qscore_cutoff(seq: str, cutoff: float): -> bool
    '''Function that takes a sequence, in this case of phred quality scores,
        and returns True if each score in the sequence meets the given cutoff, False if not'''
    For each char in seq, call bioinfo.convert_phred and compare the converted score to the given cutoff. If any of them are less than the cutoff, return False. If you get to the end of the string, return True.
Input: seq: 'AAAEEE', cutoff: 25.0
Expected Output: True

def write_record(file, record: list): 
    '''Function that takes a file header and a list of values and writes the list to the file.'''
    use with open on the given file, and use file.writelines on the list to write each element of the list to the file. 
Input: R1_Index1.fastq, ['Record1_header\n', 'Record1_seq\n', '+\n', 'Record1_qscores\n']
Expected Output: The file R1_Index1.fastq has the 4 elements of the list written on separate lines
```

create dictionary of all swapped and matching known-index pairs (576 total):
```
i.e. {'GTAGCGTA-GTAGCGTA': 0, 'GTAGCGTA-AACAGCGA': 0, ...} 
* use itertools
``` 

create dictionary of filehandles for all output files:
```
i.e. {"R1-B1": open("R1-B1.fastq", "w"), "R1-B9": open("R1-B9.fastq", "w"), ...}
```
This is to avoid opening and closing the files every time you want to add a record. That would not be feasible for 360million+ records. Instead, open them all at the beginning and keep them open throughout the loop. Need to remember to close all the files at the end.

```
1. Initialize unknown count, The 4 input files, and the qscore cutoff.

2. Open all four of the input files (R1, R2, R3, R4).

3. Use a while True loop, and read each line in the records from R1 and R4, as well as the index and qscore lines from R2 and R3. If the header line from R1 is an empty string, you know you have hit the end of the file and can break from the while loop. 

    4. Calculate and store the reverse complement of the R3 index (using reverse complement function).

    5. If 'r2index-rev_comp_r3index' not in dictionary of known index pairs, write associated reads from R1 and R4 to 'unknownR1' and 'unknownR4' output files with indexes appended to the header line (using write_lines function). Increment unknown count.

    6. Else, 
        if meets_qscore_cutoff returns False, write the associated reads from R1 and R4 to 'unknownR1' and 'unknownR4' output files with indexes appended to the header line (using write_lines function). Increment unknown count.

        7. Else, increment count of 'index from R2-reverse complement' position in the dictionary. 

        8. If reverse complement == index from R2, write associated reads from R1 and R4 to respective matching index-pair output file with indexes appended to the header line (using write_lines function). 

        9. Else, write associated reads from R1 and R4 to respective hopped index output files with indexes appended to the header line (using write_lines function). 

10. After the loop, loop through the dictionary of fileheaders and close all the files.

11. Use dictionary of index pairs to get the sums of the matched index pairs and the hopped index pairs. Loop through each key value pair, and split the key string on the '-'. If the 2 elements are equal, then it is a matched key, so increase the matched sum by the value at that key. Else, it's a hopped key, so increase the hopped sum by that value. Return the unkown, matched, and hopped sums at the end. 
```

