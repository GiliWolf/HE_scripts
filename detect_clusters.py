"""
Author: Gili Wolf
Date: 02-03-2024
Levanon Lab
----------------------------
This script is designed to analyze and detect the editing events from the overall general mismatches events.
the script output a CSV file with the following parameters:
    1. Read_ID 
    2. Chromosome: 
    3. Position:
    4. Alignment_length:
    5. Reference_Sequence:
    6. cigar:
    7. blocks:
    8. HE_Sequence:
    9. Number_of_MM:
    10. Number_of_Editing_Sites:
    11. Editing_Sites_List:


Usage: python detect_clusters.py <bam_path> <fasta_path> <output_path> <ref_base> <alt_base>
 
"""

import pysam
import csv
import sys

if len(sys.argv) != 6:
    print("Usage: python detect_clusters.py <bam_path> <fasta_path> <output_path> <ref_base> <alt_base>")
    sys.exit(1)

bam_path = sys.argv[1]
fasta_path = sys.argv[2]
output_path = sys.argv[3]
ref_base = sys.argv[4]
alt_base = sys.argv[5]

# Open the BAM file for reading
bam_file = pysam.AlignmentFile(bam_path, "rb")

# Open the FASTA genome file
fasta_file = pysam.FastaFile(fasta_path)

# Open the output file for writing
output_file = open(output_path, "w", newline='')
csv_writer = csv.writer(output_file)

# Write the header row to the CSV file
header = ['Read_ID', 'Chromosome', 'Position','Alignment_length', 'Reference_Sequence', 'cigar','blocks', 'HE_Sequence', 'Number_of_MM', 'Number_of_Editing_Sites', 'Editing_Fracture', 'Editing_Sites_List']
csv_writer.writerow(header)

# Iterate through each read in the BAM file
for read in bam_file:
    # Check if the read is mapped
    if not read.is_unmapped:
        # Get the read's mapped region coordinates 
        chromosome = bam_file.get_reference_name(read.reference_id)  # chromosome name
        position = read.reference_start  # start position of the allignment
        
        # Extract the sequence from the genome:
        # get the cigar
        cigar = read.cigarstring
        # get blocks of alligned regions (if size of blocks > 1 -> the read is spliced)
        blocks = read.get_blocks()
        reference_sequence=""
        # iterate over the blocks of refrence sequence using block[0](start position) and block[1](end position)
        for block in blocks:
            reference_sequence += str(fasta_file.fetch(chromosome, block[0], block[1])).upper()
        
        allignment_length = len(reference_sequence)

        # get the sequence of the read itself.
        # if reversed - get the reverse complemented (get_forward_sequence)
        HE_sequence = read.get_forward_sequence() if read.is_reverse else read.query_alignment_sequence

        # list of the quality of each base
        quality = read.query_alignment_qualities
        
        # MM parameters:
        num_of_mm = 0 # total MM
        num_of_editing_sites = 0 # editing sites (based on the current file base combination)
        editing_sites = {} # map : (position: quality)

        # Extract MM and editing sites:
        for i in range(min(len(reference_sequence), len(HE_sequence), len(quality))):
            temp_ref_base = reference_sequence[i] #reference base of this read
            temp_alt_base = HE_sequence[i] # alternate base of this read
            if temp_ref_base != temp_alt_base: # MM
                num_of_mm += 1
                if (temp_ref_base == ref_base) and (temp_alt_base == alt_base): # MM is same as base combination
                   num_of_editing_sites += 1
                   editing_sites[i] = quality[i] # key: position, value: position's quality

        editing_fracture = 0 if num_of_mm == 0 else num_of_editing_sites / num_of_mm

        # Write the parameters to the CSV file:
        # 'Read_ID', 'Chromosome', 'Position','Alignment_length', 'Reference_Sequence', 'cigar','blocks', 'HE_Sequence', 'Number_of_MM', 'Number_of_Editing_Sites', 'editing_Fracture', 'Editing_Sites_List'
        csv_writer.writerow([read.query_name, chromosome, position, allignment_length, reference_sequence,cigar, blocks, HE_sequence, num_of_mm, num_of_editing_sites,editing_fracture, editing_sites])

# Close all of the files
bam_file.close()
fasta_file.close()
output_file.close()