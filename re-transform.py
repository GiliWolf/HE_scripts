"""
Author: Gili Wolf
Date: 01-02-2024
Levanon Lab
----------------------------
This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) 
and insert them into the aligned sequences of a BAM file, which contained transformed sequences.

Usage: python re-transform.py <sam_file> <output_path> <original_fastq_directory>

algo: 
    1. extract reads_id and original sequences from each mate's original fastqs (in the original_reads_dir)
    2. iterate over each transformed mapped read in the sam file:
        a) check if the SAM read is mate1/mate2
        b) find matching original fastq read from (1)
        c) get original sequence and write the SAM read with the original sequence to the output file 
"""

    # quay.io/biocontainers/pysam key functions:
    #     get_query_names(self)
    #     is_paired
    #     is_read1
    #     is_read2
    #     is_reverse
    #     mate_is_reverse
    #     mate_is_unmapped
    #     is_proper_pair ?
    #     mapping_quality
    #     query_alignment_sequence / query_sequence
    #     query_name

    #     pysam.FastxFile


import pysam
import sys
import os

# arguments:
input_sam_file = sys.argv[1] 
output_sam_file = sys.argv[2]
fastq_directory = sys.argv[3]

# List all files in the directory
files = os.listdir(fastq_directory)

# Filter files ending with "mate1" and "mate2"
original_mate_1 = next((f for f in files if f.endswith("mate1")), None)
original_mate_2 = next((f for f in files if f.endswith("mate2")), None)

# Create a dictonary to store the read names (as keys) and sequences (as value) from mate1 and mate 2 Fastqs
mate1_seqs = {}
mate2_seqs = {}

# Open mate1's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(os.path.join(fastq_directory, original_mate_1)) as mate1:
    for read in mate1:
        mate1_seqs[read.name] = read.sequence

# Open mate2's unmapped Fastq original file and store all its read names in a set
with pysam.FastxFile(os.path.join(fastq_directory, original_mate_2)) as mate2:
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
