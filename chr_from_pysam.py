import pysam

# Open the BAM file
bam_file = pysam.AlignmentFile("/private10/Projects/Gili/HE_workdir/MAIN_tests/MAIN_SEARCH_and_detect/re-transform/A2C/A2C_SRR11548780_small.bam", "rb")

# Iterate through the reads
for read in bam_file:
    # Extract the chromosome (reference name)
    chromosome = read.reference_name
    print(chromosome)
    break