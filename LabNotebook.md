# Demultiplex Assignment

## Assignment the First

### Description
Separate read data of flow cell results by index pairs. If the index from read 1 and read 2 is the same, and it is a known index, it is a matched pair. If the indexes from read 1 and read 2 are both known, but they do not match, it is a hopped index pair. And if either of the indexes are not known, or if their quality is too low, they are considered unknown. Format the output into separate files, one for each known matched-index pair for each read, one for all hopped indexes for each read, and one for all unknown indexes for each read. Also determine the counts of each category of output(matched, hopped, unknown). 

### 07/21/2026

Working with Gabe

After looking at the files, we can see that the headers match in each file, indicating that they all correlate to each other. We suspect that the read fastqs have been trimmed to remove both the i5 and i7 sections from the beginning and end, as well as the indexes from the beginning and end. It looks like all the indexes start with an N, but the rest of the bases line up with a barcode (if it matched). However, when looking at the indexes from the reverse read, they are reverse complements of the indexes from the forward reads (besides the N at the beginning). 

Follow-up: Not all the indexes start with an N, actually 359 million of the 363 million do not. The strategy will be to automatically put any index with an N and associated record into the unknown output files. 

Drafted psuedocode for reverse_complement function and demultiplex logic. Decided on a dictionary for the index pairs, and a dictionary for all the filehandles. 

