from Bio import SeqIO
from Bio.Seq import Seq
import sys

# this function changed each ref_base to alt_base in all of the sequancec of a given FASTQ file 
def process_fastq(input_file, output_file, ref_base, alt_base):
    # get upper and lower case of the bases
    lower_ref_base = str(ref_base).lower()
    lower_alt_base = str(alt_base).lower()

    upper_ref_base = str(ref_base).upper()
    upper_alt_base = str(alt_base).upper()

    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        for record in SeqIO.parse(in_file, 'fastq'):
            record.seq = Seq((str(record.seq)).replace(lower_ref_base, lower_alt_base))
            record.seq = Seq((str(record.seq)).replace(upper_ref_base, upper_alt_base))
            SeqIO.write(record, out_file, 'fastq')
    
  

if __name__ == "__main__":
    input_fastq = sys.argv[1] 
    output_fastq = sys.argv[2]
    ref_base = sys.argv[3]
    alt_base = sys.argv[4]

    process_fastq(input_fastq, output_fastq, ref_base, alt_base)
