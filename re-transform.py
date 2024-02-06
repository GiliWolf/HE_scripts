"""
Author: Gili Wolf
Date: 01-02-2024
----------------------------
This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) 
and insert them into the aligned sequences of a BAM file, which contained transformed sequences.

Usage: python re-transform.py <sam_file> <output_path> <mate1> [mate2]
"""

import pysam
import sys

# arguments:
input_sam_file = sys.argv[1] 
output_sam_file = sys.argv[2]
original_mate_1 = sys.argv[3]
original_mate_2 = sys.argv[4]

# Create a dictonary to store the read names (as keys) and sequences (as value) from mate1 and mate 2 Fastqs
mate1_seqs = {}
mate2_seqs = {}

# Open mate1's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(original_mate_1) as mate1:
    for read in mate1:
        mate1_seqs[read.name] = read.sequence

# Open mate2's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(original_mate_2) as mate2:
    for read in mate2:
        mate2_seqs[read.name] = read.sequence

# change the read seq (from sam file) to the seq extracted from the same mate's fastq file, based on read's name 
def getOriginal(read, seq_dict):
    # Access the name of the read
    transformed_id = read.query_name

    # Check if the transformed_id is in the set of mate1 read names
    if transformed_id in seq_dict:
        # get original sequence
        original_seq = seq_dict[transformed_id]
        # modify the sam file sequence to the original 
        read.query_sequence = original_seq
        
        # write it to the output SAM file
        output_sam.write(read)


# Open the SAM file for reading
with pysam.AlignmentFile(input_sam_file, "rb") as samfile:
    # Open a new SAM file for writing
    with pysam.AlignmentFile(output_sam_file, "wh", header=samfile.header) as output_sam:
        # for each read in the sam file
        for read in samfile:
            # get original read of mate1
            if read.is_read1: 
                getOriginal(read, mate1_seqs)
            # get original read of mate2
            if read.is_read2: 
                getOriginal(read, mate2_seqs)
