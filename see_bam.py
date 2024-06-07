import pysam

bam_file = "/private10/Projects/Gili/HE_workdir/first_part/PE_test3/second_map/G2A/SRR11548778_G2A_Aligned.out.bam"

# Open the BAM file
with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam:
        print(read)