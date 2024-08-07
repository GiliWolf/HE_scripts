## TO-DO - 
#   Multimappers -
#       (1) Choose read with top fraction 
#       (2) add to json *uniqluey* before reads (group by ID)


"""
Author: Gili Wolf
Date: 02-03-2024
Levanon Lab
----------------------------
This script is designed to analyze and detect the editing events from the overall general mismatches events.
the script output a CSV file with the following parameters:
    1. Read_ID 
    2. Chromosome 
    3. Strand (+/-)
    4. Position (0 based)
    5. Alignment_length (clipping isn't included, length of the *coding* regions in case of splicing)
    6. Read_Sequence (original ran read)
    7. Visualize Allignment ('*': editing site, 'X': other MM, '|': match)
    8. Reference_Sequence (clipping isn't included, introns are not incloded, from the genome fasta file)
    9. cigar
    10. SAM file flags
    11. Genomic_Position_Splicing_Blocks (coding regions in case of splicing (0-based))
    12. Read_Relative_Splicing_Blocks (splicing blocks related to the read's position (0-based))
    13. Number_of_Editing_Sites (number of mismatches of the ref_base to alt_base kind)
    14. Number_of_MM (number of total mismatches of Read_Sequence compare to Reference_Sequence)
    15. EditingSites_to_PhredScore_Map (key: position, value: phred quality score)
    16. MM_to_PhredScore_Map (key: position, value: phred quality score)
----------------------------

usage: detect_clusters.py [-h] -i BAM_PATH -g FASTA_PATH -o OUTPUT_PATH -rb REF_BASE -ab ALT_BASE [-c {all,basic}]

  -h, --help            show this help message and exit
  -i BAM_PATH, --bam_path BAM_PATH
                        Path to the BAM file.
  -g FASTA_PATH, --fasta_path FASTA_PATH
                        Path to the FASTA genome file.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to the output CSV file.
  -rb REF_BASE, --ref_base REF_BASE
                        Reference base.
  -ab ALT_BASE, --alt_base ALT_BASE
                        Alternate base.
  -c {all,basic}, --output_columns {all,basic}
                        Choose the type of output columns: 'all' or 'basic' ('Read_ID', 'Chromosome',
                        'Position','Alignment_length','Read_Sequence', 'Reference_Sequence','Number_of_MM', 'Number_of_Editing_Sites')
 
"""

import pysam
import csv
import argparse
#from Bio.Seq import Seq


# Parse command line arguments
parser = argparse.ArgumentParser(description="This script is designed to analyze and detect the editing events from the overall general mismatches events.")
parser.add_argument("-i","--bam_path", type=str, required=True, help="Path to the BAM file.")
parser.add_argument("-g", "--fasta_path", type=str, required=True, help="Path to the FASTA genome file.")
parser.add_argument("-o","--output_path", type=str, required=True, help="Path to the output CSV file.")
parser.add_argument("-rb", "--ref_base", type=str, required=True, help="Reference base.")
parser.add_argument("-ab", "--alt_base", type=str, required=True, help="Alternate base.")
parser.add_argument("-c","--output_columns", type=str, choices=["all", "basic"], default='all', help="Choose the type of output columns: 'all' or 'basic' ('Read_ID', 'Chromosome', 'strand', 'Position','Alignment_length','Read_Sequence', 'Reference_Sequence','Number_of_MM', 'Number_of_Editing_Sites', 'Editing_to_Total_MM_Fraction')")

args = parser.parse_args()
# Extract command line arguments
bam_path = args.bam_path
fasta_path = args.fasta_path
output_path = args.output_path
ref_base = args.ref_base
alt_base = args.alt_base
output_columns = args.output_columns

# Open the BAM file for reading
bam_file = pysam.AlignmentFile(bam_path, "r", check_sq = False)

# Open the FASTA genome file
fasta_file = pysam.FastaFile(fasta_path)

# base combination to MM column name map
col_names_MM_map = {
    ('A', 'C'): 'A2C_MM',
    ('A', 'G'): 'A2G_MM',
    ('A', 'T'): 'A2T_MM',
    ('C', 'A'): 'C2A_MM',
    ('C', 'G'): 'C2G_MM',
    ('C', 'T'): 'C2T_MM',
    ('G', 'A'): 'G2A_MM',
    ('G', 'C'): 'G2C_MM',
    ('G', 'T'): 'G2T_MM',
    ('T', 'A'): 'T2A_MM',
    ('T', 'C'): 'T2C_MM',
    ('T', 'G'): 'T2G_MM',
}
mm_col_names = list(col_names_MM_map.values()) + ['Ref2N_MM', 'NtoAlt_MM']


# keep all the rows to be written at the end of the iteration
rows = []

def identify_Editing_Site(read_ref_base, read_alt_base):
    # EDITING SITE:
    return (read_ref_base == ref_base) and (read_alt_base == alt_base)


def identify_MM(read_ref_base, read_alt_base, quality, index, editing_sites, detected_MM_map, all_MM):
    #Unknown MM:
    # alt base is N
    if read_alt_base == 'N':
        detected_MM_map['Ref2N_MM'].append(index)
        all_MM[index] = quality[index] # add to all MM quality map
        return 'X'
    # ref base is N
    if read_ref_base == 'N':
        detected_MM_map['NtoAlt_MM'].append(index)
        all_MM[index] = quality[index] # add to all MM quality map
        return 'X'

    # EDITING SITE:
    if (read_ref_base == ref_base) and (read_alt_base == alt_base):
        editing_sites[index] = quality[index]
        detected_MM_map[col_names_MM_map[(read_ref_base, read_alt_base)]].append(index)
        return '*'
    
    # known MM
    detected_MM_map[col_names_MM_map[(read_ref_base, read_alt_base)]].append(index)
    all_MM[index] = quality[index] # add to all MM quality map
    return 'X'

def get_read_blocks(genomic_blocks):
    read_blocks = []
    offset = 0
    for block in genomic_blocks:
        end = offset + (block[1] - block[0]) - 1
        read_blocks.append((offset, end))
        offset = end + 1
    return read_blocks
    

# Iterate through each read in the BAM file
for read in bam_file:
        # keep track of all MM detected in the read
        detected_MM_map = {key: [] for key in mm_col_names}
        # Check if the read is mapped
        if not read.is_unmapped:
            # BASIC INFO:
            # Get the read's mapped region coordinates 
            chromosome = bam_file.get_reference_name(read.reference_id)  # chromosome name
            position = read.reference_start
            genomic_blocks = read.get_blocks()
            
            # iterate over the genomic_blocks to:
            #   1. extract refrence sequence using block[0](start position) and block[1](end position)
            #   2. extract read_relative_splicing_blocks  - how the read is being splices based on read's positions

            reference_sequence = "".join([fasta_file.fetch(chromosome, block[0], block[1]).upper() for block in genomic_blocks])
            
            allignment_length = len(reference_sequence)
            # if reversed - get the reverse complemented of the alignment
            complementary = { 'A':'T', 'T':'A', 'G':'C','C':'G', 'N':'N'}
            rev_comp_seq = ''.join(([complementary[i] for i in reversed(reference_sequence)]))
            #reference_sequence = str(rev_comp_seq) if read.is_reverse else reference_sequence
            # get the sequence of the read itself.
            read_sequence = read.query_alignment_sequence
            
            # get strand orientation - +: forward, -: reverse
            strand = '-' if read.is_reverse else '+'

            # initilize MM parameters:
            num_of_mm = 0 # total MM
            num_of_editing_sites = 0 # editing sites (based on the current file base combination)
            editing_sites = {} # map : (position: quality)

            
            # ADDITIONAL
            quality = read.query_alignment_qualities
            flag = read.flag
        
            editing_sites = {} # map : (position: quality)

            # visualization of the allignment - '*': editing site, 'X': other MM, '|': match.
            visualize_allignment = ''

            # init all MM map 
            all_MM = {}
            # Extract MM and editing sites:
            for i in range(min(len(reference_sequence), len(read_sequence), len(quality))):
                read_ref_base = reference_sequence[i] #reference base of this read
                read_alt_base = read_sequence[i] # alternate base of this read
                if read_ref_base != read_alt_base: # MM
                    num_of_mm += 1
                    visualize_allignment += identify_MM(read_ref_base, read_alt_base, quality, i, editing_sites,detected_MM_map, all_MM)
                    num_of_editing_sites = len(editing_sites)
                else:
                    visualize_allignment += '|'
                
            editing_fraction = (num_of_editing_sites/num_of_mm) if num_of_mm != 0 else 1
            # append the parameters to the CSV file rows:
            row = [
                read.query_name, chromosome, position, flag, read_sequence,visualize_allignment, reference_sequence, num_of_editing_sites, num_of_mm, editing_fraction
            ]

            rows.append(row)



with open(output_path, 'w', newline='') as output_file:
    csv_writer = csv.writer(output_file, csv.QUOTE_MINIMAL)
    header = ['Read_ID','chromosome', 'position', 'flag', 'Read_Sequence','visualize_allignment' 'Reference_Sequence','Number_of_Editing_Sites','Number_of_MM', 'editing_fraction']

        
    csv_writer.writerow(header)
    csv_writer.writerows(rows)

# Close all of the files
bam_file.close()
fasta_file.close()
