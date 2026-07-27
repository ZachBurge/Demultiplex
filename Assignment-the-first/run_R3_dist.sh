#!/bin/bash

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1

R3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

/usr/bin/time -v pixi run python create_qscore_distribution.py -i $R3 -o R3 -l 8