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
        a) check if the SAM read is mate1/mate2
        b) find matching original fastq read from (1)
        c) get original sequence and write the SAM read with the original sequence to the output file 


usage: re-transform.py [-h] -s INPUT_SAM_FILE -o OUTPUT_SAM_FILE -i FASTQ_PATH -I FASTQ_PATH_MATE2 -pe {0,1}

Required Arguments:
  -h, --help                show this help message and exit
  -s --input_sam_file       Input SAM file.
  -o --output_sam_file      Output SAM file.
  -i --fastq_path           Path to FASTQ file (mate1).
  -I --fastq_path_mate_2    Path to FASTQ file for mate2.
  -pe --pair_end {0,1}      Specify if single-end (0) or paired-end (1) data.
"""

import pysam
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to recover the original sequences from one or two Fastq files (for paired-end sequencing) and insert them into the aligned sequences of a BAM file, which contains transformed sequences.")

parser.add_argument("-s", "--input_sam_file",dest ="input_sam_file", type=str, help="Input SAM file.", required=True)
parser.add_argument("-o","--output_sam_file",dest="output_sam_file", type=str, help="Output SAM file.",required=True)
parser.add_argument("-i", "--fastq_path", dest="fastq_path", type=str, help="Path to FASTQ file.",required=True)
parser.add_argument("-I", "--fastq_path_mate_2", dest="fastq_path_mate2", type=str, help="Path to FASTQ file for mate2.",required=False)
parser.add_argument("-pe","--pair_end",dest = "pair_end", type=int, choices=[0, 1], help="Specify if single-end (0) or paired-end (1) data.",required=True)

args = parser.parse_args()

# get args
input_sam_file = args.input_sam_file
output_sam_file = args.output_sam_file
mate1_path = args.fastq_path
pair_end = args.pair_end
if pair_end:
    if args.fastq_path_mate2:
        mate2_path=args.fastq_path_mate2
    else:
        print("re-transform.py: the following arguments is required for PE data: -I/--fastq_path_mate_2")
        exit()



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
        original_seq = seq_dict[transformed_id]
        # modify the sam file sequence to the original 
        read.query_sequence = original_seq
        # restore the read qualities
        read.query_qualities = qualities[transformed_id]
        # write it to the output SAM file
        output_sam.write(read)


# Open the SAM file for reading
try:
    with pysam.AlignmentFile(input_sam_file, "rb") as samfile:
        # Open a new SAM file for writing
        with pysam.AlignmentFile(output_sam_file, "wh", header=samfile.header) as output_sam:
            # for each read in the sam file
            for read in samfile:
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
except Exception as e:
    print(f"RE-TRANSFORM: error parsing files {e}")
    exit()

def benchmark_write(samfile, threads):
    bamfile = samfile[:-4]+'_pysam.bam'
    infile = pysam.AlignmentFile(samfile, 'r', threads=threads)
    outfile = pysam.AlignmentFile(bamfile, mode='wb',
                                    template=infile, threads=threads)
    for s in infile:
        outfile.write(s)
        
def process_reads_batch(reads, fasta_file, bam_file):
    rows = []
    for read in reads:
        row = process_read(read, fasta_file, bam_file)
        rows.append(row)
    return rows

def main():
    # Open the BAM file for reading
    bam_file = pysam.AlignmentFile(bam_path, "rb")

    # Open the FASTA genome file
    fasta_file = pysam.FastaFile(fasta_path)

    # Keep all the rows to be written at the end of the iteration
    rows = []

    # Read batches of reads
    batches = []
    batch = []

    # split file into batches of reads
    for read in bam_file:
        batch.append(read)
        if len(batch) == batch_size:
            batches.append(batch)
            batch = []
    if batch:
        batches.append(batch)  # Add the remaining reads

    # Use ThreadPoolExecutor to parallelize processing
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [executor.submit(process_reads_batch, batch, fasta_file, bam_file) for batch in batches]
        for future in as_completed(futures):
            rows.extend(future.result())

    # Write the output to the CSV file
    with open(output_path, 'w', newline='') as output_file:
        csv_writer = csv.writer(output_file, csv.QUOTE_MINIMAL)
        if output_columns == "basic": # basic columns
            header = ['Read_ID', 'Chromosome', 'Strand', 'Position_0based', 'Alignment_length', 'Read_Sequence', 'Reference_Sequence', 'Number_of_MM', 'Number_of_Editing_Sites']
        else: # all columns
            header = ['Read_ID', 'Chromosome', 'Strand', 'Position_0based', 'Alignment_length', 'Read_Sequence', 'Visualize_Allignment', 'Reference_Sequence', 'cigar', 'flag', 'Genomic_Position_Splicing_Blocks_0based', 'Read_Relative_Splicing_Blocks_0based', 'Number_of_Editing_Sites', 'Number_of_total_MM', 'EditingSites_to_PhredScore_Map', 'MM_to_PhredScore_Map']
            header.extend(mm_col_names)
            
        csv_writer.writerow(header)
        csv_writer.writerows(rows)

    # Close all of the files
    bam_file.close()
    fasta_file.close()

if __name__ == "__main__":
    main()