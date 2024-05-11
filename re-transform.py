"""
Author: Gili Wolf
Date: 01-02-2024
Levanon Lab
----------------------------
This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) 
and insert them into the aligned sequences of a BAM file, which contained transformed sequences.

Usage: python re-transform.py <sam_file> <output_path> <original_fastq_path> <pair_end[0/1]>

algo: 
    1. extract reads_id and original sequences from each mate's original fastqs (in the original_reads_dir)
    2. iterate over each transformed mapped read in the sam file:
        a) check if the SAM read is mate1/mate2
        b) find matching original fastq read from (1)
        c) get original sequence and write the SAM read with the original sequence to the output file 
"""

import pysam
import sys
import os
import argparse

# Controling Arguments:
parser = argparse.ArgumentParser(description="Usage: python re-transform.py -s <sam_file> -o <output_path> -i <original_fastq_path> [-I <original_fastq_path_2>] <pair_end[0/1]>")

parser.add_argument("-s", "--input_sam_file",dest ="input_sam_file", type=str, help="Input SAM file.")
parser.add_argument("-o","--output_sam_file",dest="output_sam_file", type=str, help="Output SAM file.")
parser.add_argument("-i", "--fastq_path", dest="fastq_path", type=str, help="Path to FASTQ file.")
parser.add_argument("-I", "--fastq_path_mate_2", dest="fastq_path_mate2", type=str, help="Path to FASTQ file for mate2.")
parser.add_argument("-pe","--pair_end",dest = "pair_end", type=int, choices=[0, 1], help="Specify if single-end (0) or paired-end (1) data.")

args = parser.parse_args()

# get args
input_sam_file = args.input_sam_file
output_sam_file = args.output_sam_file
mate1_path = args.fastq_path
pair_end = args.pair_end
if pair_end:
    mate2_path=args.fastq_path_mate2



# Create a dictonary to store the read names (as keys) and sequences (as value) from mate1 and mate 2 Fastqs
mate1_seqs = {}
mate2_seqs = {}

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
        original_seq = seq_dict[transformed_id]
        # modify the sam file sequence to the original 
        read.query_sequence = original_seq
        #restore the read qualities
        read.query_qualities = qualities[transformed_id]
        # write it to the output SAM file
        output_sam.write(read)


# Open the SAM file for reading
with pysam.AlignmentFile(input_sam_file, "rb") as samfile:
    # Open a new SAM file for writing
    with pysam.AlignmentFile(output_sam_file, "wh", header=samfile.header) as output_sam:
        # for each read in the sam file
        for read in samfile:
            #bug in pysam- doesn't pass the read qualities when passed to a function
            # solution: pass the quailities seperatly and restore them in the getOriginal function
            qualities = read.query_qualities
            if pair_end:
                # get original read of mate1
                if read.is_read1: 
                    getOriginal(read, mate1_seqs, mate1_qualities)
                # get original read of mate2 (if pair_end == 1)
                if read.is_read2: 
                    getOriginal(read, mate2_seqs, mate2_qualities)
            else:
                getOriginal(read, mate1_seqs, qualities)

