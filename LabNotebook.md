# Demultiplex Assignment

## Assignment the First

### Description
Separate read data of flow cell results by index pairs. If the index from read 1 and read 2 is the same, and it is a known index, it is a matched pair. If the indexes from read 1 and read 2 are both known, but they do not match, it is a hopped index pair. And if either of the indexes are not known, or if their quality is too low, they are considered unknown. Format the output into separate files, one for each known matched-index pair for each read, one for all hopped indexes for each read, and one for all unknown indexes for each read. Also determine the counts of each category of output(matched, hopped, unknown). 

### 07/21/2026

Working with Gabe

After looking at the files, we can see that the headers match in each file, indicating that they all correlate to each other. We suspect that the read fastqs have been trimmed to remove both the i5 and i7 sections from the beginning and end, as well as the indexes from the beginning and end. It looks like all the indexes start with an N, but the rest of the bases line up with a barcode (if it matched). However, when looking at the indexes from the reverse read, they are reverse complements of the indexes from the forward reads (besides the N at the beginning). 

Follow-up: Not all the indexes start with an N, actually 359 million of the 363 million do not. The strategy will be to automatically put any index with an N and associated record into the unknown output files. 

Drafted psuedocode for reverse_complement function and demultiplex logic. Decided on a dictionary for the index pairs, and a dictionary for all the filehandles. 

### 07/25/2026

### Part 1
Initial data analysis

File R1 and R4 contain the reads, and R2 and R3 contain the indexes. 

| File | Data | Length of reads | Number of records | Number of unique barcodes | Number of barcodes with "N's" |
| ---- | ---- | ---- | ---- | ---- | ---- |
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | 363246735 | | |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | 363246735 | 69699 | 3976613 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | 363246735 | 53613 | 3280930 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | 363246735 | | |

Command used to get length of reads:
```
zcat 1294_S1_L008_R[1234]_001.fastq.gz | head -4 | grep -A1 "^@" | tail -1 | wc -c
```
Need to subtract one from the result because it is including the newline character at the end.

Command used to get number of records: 
```
zcat 1294_S1_L008_R[1234]_001.fastq.gz | grep "^@" | wc -l
```

Command used to get number of unique barcodes:
```
zcat 1294_S1_L008_R[23]_001.fastq.gz | grep -A1 "^@" | grep -v "^@" | grep -v "^--" | sort | uniq -c | wc -l
```

Command used to get number of barcodes with N's:
```
zcat 1294_S1_L008_R[23]_001.fastq.gz | grep -A1 "^@" | grep -v "^@" | grep -v "^--" | grep "N" | wc -l
```

The phred encoding for these data is Phred+33 because there are hashtags in the quality score strings, which are only used in phred+33.

### 07/26/2026

### Part 1

Realized that my command was wrong for determining number of barcodes that had an N. I was counting only the ones that started with an N. Interestingly, when I fixed it, there were only about 200 more in R2, meaning most of the barcodes that have an N start with an N.

Created script for making the qscore distributions. My original plan was to make a 2D numpy array with 101 rows and 363,246,735	columns in each row, but after building the script and testing it, numpy told me that was too big. An array that size would take 250+ Gbs. I still wanted to use numpy arrays for this, so I decided to create a 2D array to count the frequencies of each quality score at each nucleotide position. So there are still 101 rows, but now there are only 42 columns, because the maximum observed quality score in the files is 'J', which equates to a quality score of 41. Now when looping through the file, everytime a quality score is observed, the count at that position in the array is incremented. 

To calculate the mean of each nucleotide position, I created a second numpy array that was just a range from 0 to 41. I then multiplied each row in the 2D array by this array, and summed each total. This gave me an array of length 101, with each element being the sum of the observed quality scores across all the reads at that nucleotide position. All I had to do after that was divide this array by the number of records in the file, which is a known value (363,246,735) to get an array of all the mean quality scores at each nucleotide position. 

Since these are very large files, I created 4 bash scripts to run simultaneously as sbatch jobs. I also used argparse to make the script dynamic enough to work for all 4 files separately. 

Here are the /usr/bin/time results from the sbatch jobs:
| Job | Wall Time | CPU | MRSS | Exit Status |
| ---- | ---- | ---- | ---- | ---- |
| run_R1_dist.sh | 2:14:42 | 99% | 81.052 MB | 0 |
| run_R2_dist.sh | 15:05.23 | 99% | 79.636 MB | 0 |
| run_R3_dist.sh | 15:09.41 | 99% | 79.332 MB | 0 |
| run_R4_dist.sh | 3:30:53 | 99% | 83.636 MB | 0 | 