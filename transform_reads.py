from Bio import SeqIO
from Bio.Seq import Seq

def replace_a_with_g(sequence):
    """Replace all 'A' with 'G' in a given sequence."""
    return sequence.replace('A', 'G')

def process_fastq(input_file, output_file):
    """Reads a FASTQ file, replaces 'A' with 'G' in each read, and writes to a new FASTQ file."""
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        i=0
        for record in SeqIO.parse(in_file, 'fastq'):
            i=i+1
            record.seq = Seq((str(record.seq)).replace('A', 'G'))
            SeqIO.write(record, out_file, 'fastq')
            if i==100:
                break
  

if __name__ == "__main__":
    input_fastq = "/private10/Projects/Gili/HE_workdir/first_part/SRA2/SRR11548778_1.fastq"  # Replace with your input FASTQ file
    output_fastq = "SRR11548778_1_transformed.fastq"  # Replace with the desired output FASTQ file

    process_fastq(input_fastq, output_fastq)
