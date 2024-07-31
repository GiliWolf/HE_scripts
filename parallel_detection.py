import csv
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import argparse
import threading
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
Parallel Processing:

The BAM file will be divided into batches for parallel processing. 
The size of each batch is determined as follows:
    (-) If the batch size (-b) is provided, it will always be used.
    (-) Otherwise, if the number of reads (-n) is provided, each batch will be sized at int(num_of_reads / max_threads).
    (-) If neither the batch size nor the number of reads is provided, a default batch size of 50 will be used.
----------------------------

usage: detect_clusters.py [-h] -i BAM_PATH -g FASTA_PATH -o OUTPUT_PATH -rb REF_BASE -ab ALT_BASE [-c {all,basic}]

  -h --help            show this help message and exit
  -i --bam_path 
                        Path to the BAM file.
  -f --fasta_path 
                        Path to the FASTA genome file.
  -I --bam_index_path 
                        Path to the index BAI file.
  -F --fasta_index_path 
                        Path to the index FAI file.
  -o --output_path 
                        Path to the output CSV file.
  -rb --ref_base 
                        Reference base.
  -ab --alt_base 
                        Alternate base.
  -t --threads  
                        Number of threads to parallel processing.(default: 5)
  -n --num_of_reads 
                       number of reads in the bam file. (batch size will be calculated using int(reads_count/threads))
  -c {all,basic}, --output_columns {all,basic}
                        Choose the type of output columns: 'all' or 'basic' ('Read_ID', 'Chromosome',
                        'Position','Alignment_length','Read_Sequence', 'Reference_Sequence','Number_of_MM', 'Number_of_Editing_Sites')
 
"""
# init lock for mutexing the reading of the FASTA file 
fasta_lock = threading.Lock()

# Parse command line arguments
parser = argparse.ArgumentParser(description="This script is designed to analyze and detect the editing events from the overall general mismatches events.")
parser.add_argument("-i","--bam_path", type=str, required=True, help="Path to the BAM file.")
parser.add_argument("-f", "--fasta_path", type=str, required=True, help="Path to the FASTA genome file.")
parser.add_argument("-I","--bam_index_path", type=str, required=True, help="Path to the index BAI file.")
parser.add_argument("-F","--fasta_index_path", type=str, required=True, help="Path to the index FAI file.")
parser.add_argument("-o","--output_path", type=str, required=True, help="Path to the output CSV file.")
parser.add_argument("-rb", "--ref_base", type=str, required=True, help="Reference base.")
parser.add_argument("-ab", "--alt_base", type=str, required=True, help="Alternate base.")
parser.add_argument("-t", "--threads", type=int, required=False, default=5, help="Number of threads to parallel processing.(default: 5)")
parser.add_argument("-n", "--num_of_reads", type=int, required=False, help="Specify the total number of reads in the BAM file, in order to the batch size to be calculated as int(num_of_reads/threads). If a fixed batch size is preferred, use the -b option instead.")
parser.add_argument("-b", "--batch", type=int, required=False, default=50, help="Specify the batch size, i.e., the number of reads to be processed by each thread. If the batch size should be determined by the read count, use the -n option instead (default: 50)")
parser.add_argument("-c","--output_columns", type=str, choices=["all", "basic"], default='all', help="Choose the type of output columns: 'all' or 'basic' ('Read_ID', 'Chromosome', 'strand', 'Position','Alignment_length','Read_Sequence', 'Reference_Sequence','Number_of_MM', 'Number_of_Editing_Sites', 'Editing_to_Total_MM_Fraction')")

# Extract command line arguments
args = parser.parse_args()
bam_path = args.bam_path
fasta_path = args.fasta_path
bam_index_path = args.bam_index_path
fasta_index_path = args.fasta_index_path
output_path = args.output_path
ref_base = args.ref_base
alt_base = args.alt_base
output_columns = args.output_columns
max_threads = args.threads
reads_count = args.num_of_reads
user_batch_size = args.batch
default_batch_size = 50

# calculate batch size
if (user_batch_size and (user_batch_size>0)):
    batch_size = user_batch_size
elif (reads_count and (reads_count>0)):
    batch_size = int(reads_count/max_threads)
else:
    batch_size =default_batch_size

# Base combination to MM column name map
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
"""
function to identify the type of the MM: editing / base to base / unknown (N) to base / base to unknown (N)
"""
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

"""
parse the genomic blocks (genomic positions of the exons) to be refrenced as the read positions (in range of [0, seq's length])
"""
def get_read_blocks(genomic_blocks):
    read_blocks = []
    offset = 0
    for block in genomic_blocks:
        end = offset + (block[1] - block[0]) - 1
        read_blocks.append((offset, end))
        offset = end + 1
    return read_blocks

"""
process read to extract data for the output file
"""
def process_read(read, fasta_file, bam_file):
    try:
    # keep track of all MM detected in the read
        detected_MM_map = {key: [] for key in mm_col_names}
        # Check if the read is mapped
        if not read.is_unmapped:
            # BASIC INFO:
            # Get the read's mapped region coordinates 
            chromosome = bam_file.get_reference_name(read.reference_id)  # chromosome name
            position = read.reference_start  # start position of the allignment
            
            # Extract the sequence from the genome:
            # get Genomic_Position_Splicing_Blocks of alligned regions (if size of Genomic_Blocks > 1 -> the read is spliced)
            genomic_blocks = read.get_blocks()
            
            # iterate over the genomic_blocks to:
            #   1. extract refrence sequence using block[0](start position) and block[1](end position)
            #   2. extract read_relative_splicing_blocks  - how the read is being splices based on read's positions
            read_blocks = []
            with fasta_lock:
                reference_sequence = "".join([fasta_file.fetch(chromosome, block[0], block[1]).upper() for block in genomic_blocks])
            read_blocks = get_read_blocks(genomic_blocks)
            
            allignment_length = len(reference_sequence)
            # if reversed - get the reverse complemented of the alignment
            complementary = { 'A':'T', 'T':'A', 'G':'C','C':'G', 'N':'N'}
            rev_comp_seq = ''.join(([complementary[i] for i in reversed(reference_sequence)]))
            reference_sequence = str(rev_comp_seq) if read.is_reverse else reference_sequence
            # get the sequence of the read itself.
            read_sequence = read.query_alignment_sequence
            
            # get strand orientation - +: forward, -: reverse
            strand = '-' if read.is_reverse else '+'

            # initilize MM parameters:
            num_of_mm = 0 # total MM
            num_of_editing_sites = 0 # editing sites (based on the current file base combination)
            editing_sites = {} # map : (position: quality)

            
            if output_columns == "basic":
                # Extract number of MM and editing sites:
                for i in range(min(len(reference_sequence), len(read_sequence), len(quality))):
                    read_ref_base = reference_sequence[i] #reference base of this read
                    read_alt_base = read_sequence[i] # alternate base of this read
                    if read_ref_base != read_alt_base: # MM
                        num_of_mm += 1
                        # identify editing site
                        if (read_ref_base == ref_base) and (read_alt_base == alt_base):
                            num_of_editing_sites+=1

                row = [
                    read.query_name, chromosome, strand, position, allignment_length, read_sequence, reference_sequence,
                    num_of_mm, num_of_editing_sites
                ]
            else:
                # ADDITIONAL
                cigar = read.cigarstring
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

                # append the parameters to the CSV file rows:
                row = [
                    read.query_name, chromosome, strand, position, allignment_length, read_sequence,visualize_allignment, reference_sequence, cigar, flag, genomic_blocks, read_blocks, num_of_editing_sites, num_of_mm, editing_sites, all_MM
                ]
                # Convert the detected_MM_map dictionary into a list of its values and append it to the row
                row.extend([detected_MM_map[col_name] for col_name in mm_col_names])

            return row
    except Exception as e:
        print(f"error in Detection: {e}")
        
"""
process each batch of reads
"""
def process_reads_batch(reads, fasta_file, bam_file):
    rows = []
    for read in reads:
        row = process_read(read, fasta_file, bam_file)
        rows.append(row)
    return rows

def main():
    # ice cream / time - checl time
    # Open the BAM file for reading
    bam_file = pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index_path)

    # Open the FASTA genome file
    fasta_file = pysam.FastaFile(fasta_path, filepath_index=fasta_index_path)

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