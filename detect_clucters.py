import pysam
import sys
import os

bam_path = "/private10/Projects/Gili/HE_workdir/first_part/name_test/second_map/A2C/SRR11548778_A2C_Aligned.out.bam"
fasta_path = "/private10/Projects/Gili/HE_workdir/genome_setup/hg38.fa"

# Open the BAM file for reading
bam_file = pysam.AlignmentFile(bam_path, "rb")
# pysam.index(bam_file)

# Open the FASTA genome file
fasta_file = pysam.FastaFile(fasta_path)

# Iterate through each read in the BAM file
for read in bam_file:
    # Check if the read is mapped
    if not read.is_unmapped:
        # Get the read's mapped region coordinates 
        chromosome = bam_file.get_reference_name(read.reference_id)  # chromosome name
        position = read.reference_start  # coordinates
        
        # Extract the sequence from the genome using the mapped position
        reference_sequence = fasta_file.fetch(chromosome, position, position + read.query_length)
        read_sequence = read.query_sequence
        # Print or use the extracted sequence
        print(f"Read ID: {read.query_name}, Chromosome: {chromosome}, Position: {position}, reference_Sequence: {reference_sequence}, mapped_sequence = {read_sequence}")

# Close the BAM and FASTA files
bam_file.close()
fasta_file.close()
