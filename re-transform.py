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

# Check if correct number of command-line arguments are provided
if len(sys.argv) != 5:
    print("Usage: python re-transform.py <sam_file> <output_path> <original_fastq_path> <pair_end[0/1]>")
    sys.exit(1)

# arguments:
input_sam_file = sys.argv[1] 
output_sam_file = sys.argv[2]
fastq_path = sys.argv[3]
pair_end = int(sys.argv[4]) # SE: 0; PE:1

if pair_end != 0 and pair_end != 1:
    print("pair end flag value should be 0/1")
    sys.exit(1)

# # Create index file for the input SAM/BAM file
# input_index_file = input_sam_file + ".bai"
# if not os.path.exists(input_index_file):
#     pysam.index(input_sam_file)

# List all files in the directory
if (pair_end): # PE - fastq_path is a directory with both mates 
    files = os.listdir(fastq_path)
    # extract files - Filter files ending with "mate1" and "mate2"
    original_mate_1 = next((f for f in files if f.endswith("mate1")), None)
    mate1_path = os.path.join(fastq_path, original_mate_1)
    original_mate_2 = next((f for f in files if f.endswith("mate2")), None)
    mate2_path = os.path.join(fastq_path, original_mate_2)
else: #SE - fastq_path is the only mate1 file
    mate1_path = os.path(fastq_path)

# Create a dictonary to store the read names (as keys) and sequences (as value) from mate1 and mate 2 Fastqs
mate1_seqs = {}
mate2_seqs = {}

# Open mate1's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(mate1_path) as mate1:
    for read in mate1:
        mate1_seqs[read.name] = read.sequence

if (pair_end):
    # Open mate2's unmapped Fastq original file and store all its read names in a set
    with pysam.FastxFile(mate2_path) as mate2:
        for read in mate2:
            mate2_seqs[read.name] = read.sequence

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
        read.query_qualities = qualities
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
                    getOriginal(read, mate1_seqs, qualities)
                # get original read of mate2 (if pair_end == 1)
                if read.is_read2: 
                    getOriginal(read, mate2_seqs, qualities)
            else:
                getOriginal(read, mate1_seqs, qualities)

