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

After the results, the barcodes and sequences both have the lowest quality at the beginning. This aligns with what I saw previously in the initial data analysis when looking at the number of barcodes that start with an N vs the ones that just have an N anywhere - most of them that had an N it was in the first position. I think a good quality score cutoff for the indexes would be 30, since this still covers the average quality score in the first 2 positions, but it will get rid of reads that began with an N. I think it should also be 30 for the sequences for similar reasons. The worst quality scores are at the beginning of the sequence, and 30 still incorporates the average seen at those positions while cutting off anything with a higher error rate. 

### 08/03/2026

### Part 2

Finished the first draft of my demultiplex code. Here are a few ideas I had to slightly improve the efficiency of the code: 

1) Instead of calling bioinfo.convert_phred for every record in the loop, I do the calculation manually in a local function to save on function call time
2) For the hopped and unknown output files, I store the dictionary lookups in variables outside the loop so that I don't have to look them up each time. I could do this for the other 48 output files, but not sure the savings are really worth that
3) For my reverse complement function, I use the str.maketrans method. Instead of recreating the maketrans table every time the function is called, I store it in a global variable to reference within the function

Ran on my test files, and output matches what I expected!

Run with cutoff 20, took 45 mins. Distribution created looks good at first glance, there are values for each possible index pair. Matched index pairs have significantly higher counts, which is a good sign. Leslie mentioned creating a heatmap, looking into that next.

Decided on using seaborn, which is a package built on top of Matplotlib specifically for data visualization. Looking at the documentation, the implementation of a heatmap with seaborn is super simple. I created a 2D numpy array that is 24x24 (number of known indexes) to act as the ticks on the x and y axes. I also created a new dictionary called indexes_to_heatmap using the known indexes where the key was the index and the value was its position on the heatmap (0-23). Next I had to "separate" the dictionary keys (the index pairs). They are stored as "index1-index2", with the value of that key being the count of records with those index. In order to properly create the heatmap, I loop through the items of the dictionary, split the key on the '-' and store the two indexes, set the row and column values for that pair using the indexes_to_heatmap dictionary, and set the value of the numpy array at that coordinate set to be the count from the original loop. I then plot with sns.heatmap. The first time I created the heatmap, I saw a diagonal line indicating all of the matched index pairs, and everything else was dark. I realized this was because the counts of the matched index pairs was so much higher than any of the hopped index pairs. After looking into fixing this, I found matplotlib.colors LogNorm, which normalizes the data to a 0-1 scale before plotting, and after using that my heatmap looked much better. I also looked into different color schemes, and picked 'viridis'. 

### 08/04/2026

### Part 2

Decided it would be beneficial to run it with multiple different cutoffs to compare the results. Made the results and output files dynamically named, so that they would be stored separately in the /scratch folder and so that I could view them side by side. Ran it with 20, 25, 30, and no cutoff, each on their own sbatch script so they ran simultaneously.

Outputs all look somewhat similar, so adding percentage values to the distributions to more easily see the differences between them. 

### 08/06/2026

### Part 2

Leslie told me that I need to compress my output files on scratch because I am single handedly taking up 1 Tb of the 20 available Tbs for the class, since I ran it 4 times with different cutoffs. She told me to look into pigz for multithreaded zipping. Wrote some sbatch scripts for compressing the outputs with pigz, which was fairly straightforward. 
