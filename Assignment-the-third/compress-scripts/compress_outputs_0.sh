#!/bin/bash

#SBATCH --constraint=turin                # use new nodes
#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1

/usr/bin/time -v pigz -p 8 /scratch/bgmp/zburge/demux/cutoff-0/*.fastq