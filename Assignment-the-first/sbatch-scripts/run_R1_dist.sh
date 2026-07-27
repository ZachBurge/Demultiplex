#!/bin/bash

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1

R1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz

/usr/bin/time -v pixi run python create_qscore_distribution.py -i $R1 -o R1 -l 101