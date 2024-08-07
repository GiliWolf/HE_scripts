"""
Author: Gili Wolf
Date: 01-02-2024
Levanon Lab
----------------------------
This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) 
and insert them into the aligned sequences of a BAM file, which contains transformed sequences.


algo: 
    1. extract reads_id and original sequences from each mate's original fastqs (in the original_reads_dir)
    2. iterate over each transformed mapped read in the sam file:
        a) check if the BAM read is mate1/mate2
        b) find matching original fastq read from (1)
        c) get original sequence and write the BAM read with the original sequence to the output file 


usage: re-transform.py [-h] -s INPUT_BAM_FILE -o OUTPUT_BAM_FILE -i FASTQ_PATH -I FASTQ_PATH_MATE2 -pe {0,1}

Required Arguments:
  -h, --help                show this help message and exit
  -s --input_bam_file       Input BAM file.
  -o --output_bam_file      Output BAM file.
  -i --fastq_path           Path to FASTQ file (mate1).
  -I --fastq_path_mate_2    Path to FASTQ file for mate2.
  -pe --pair_end {0,1}      Specify if single-end (0) or paired-end (1) data.
"""

import pysam
import sys
import os
import argparse

# Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) and insert them into the aligned sequences of a BAM file, which contains transformed sequences.")

parser.add_argument("-s", "--input_bam_file",dest ="input_bam_file", type=str, help="Input BAM file.", required=True)
parser.add_argument("-o","--output_bam_file",dest="output_bam_file", type=str, help="Output BAM file.",required=True)
parser.add_argument("-i", "--fastq_path", dest="fastq_path", type=str, help="Path to FASTQ file.",required=True)
parser.add_argument("-I", "--fastq_path_mate_2", dest="fastq_path_mate2", type=str, help="Path to FASTQ file for mate2.",required=False)
parser.add_argument("-pe","--pair_end",dest = "pair_end", type=int, choices=[0, 1], help="Specify if single-end (0) or paired-end (1) data.",required=True)

args = parser.parse_args()

# get args
input_bam_file = args.input_bam_file
output_bam_file = args.output_bam_file
mate1_path = args.fastq_path
pair_end = args.pair_end
mate2_path = args.fastq_path_mate2 if args.fastq_path_mate2 else None

if pair_end and not mate2_path:
    print("re-transform.py: the following arguments are required for PE data: -I/--fastq_path_mate_2")
    sys.exit(1)



# Create a dictonary to store the read names (as keys) and sequences (as value) from mate1 and mate 2 Fastqs
mate1_seqs = {}
mate2_seqs = {}

# Create a dictonary to store the read names (as keys) and qualities (as value) from mate1 and mate 2 Fastqs
mate1_qualities = {}
mate2_qualities = {}



# Open mate1's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(mate1_path) as mate1:
    for read in mate1:
        mate1_seqs[read.name] = read.sequence
        mate1_qualities[read.name] = [ord(q) - 33 for q in read.quality]

if (pair_end):
    # Open mate2's unmapped Fastq original file and store all its read names in a set
    with pysam.FastxFile(mate2_path) as mate2:
        for read in mate2:
            mate2_seqs[read.name] = read.sequence
            mate2_qualities[read.name] = [ord(q) - 33 for q in read.quality]

# change the read seq (from sam file) to the seq extracted from the same mate's fastq file, based on read's name 
def getOriginal(read, seq_dict, qualities):
    # Access the name of the read
    transformed_id = read.query_name
    # Check if the transformed_id is in the set of mate1 read names
    if transformed_id in seq_dict:
        # get original sequence
        # modify the sam file sequence to the original 
        read.query_sequence = seq_dict[transformed_id]
        # restore the read qualities
        read.query_qualities = qualities[transformed_id]
        # write it to the output BAM file
        output_bam.write(read)


# Open the BAM file for reading
try:
    with pysam.AlignmentFile(input_bam_file, "rb") as bamfile:
        # Open a new BAM file for writing
        with pysam.AlignmentFile(output_bam_file, "wb", header=bamfile.header) as output_bam:
            # for each read in the sam file
            for read in bamfile:
                # bug in pysam- doesn't pass the read qualities when passed to a function
                # solution: pass the quailities seperatly and restore them in the getOriginal function
                if pair_end:
                    # get original read of mate1
                    if read.is_read1: 
                        getOriginal(read, mate1_seqs, mate1_qualities)
                    # get original read of mate2 (if pair_end == 1)
                    if read.is_read2: 
                        getOriginal(read, mate2_seqs, mate2_qualities)
                else:
                    getOriginal(read, mate1_seqs, mate1_qualities)
except:
    print("error opening files")
    exit()
