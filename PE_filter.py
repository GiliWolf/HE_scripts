from Bio import SeqIO
import pandas as pd
fastq_input_file = "/private10/Projects/Gili/HE_workdir/first_part/PE_test/first_map/SRR11548778.Unmapped.out.mate1"
he_filterd_input_file = "/private10/Projects/Gili/HE_workdir/detection/PE_output/filtered_clusters/C2A/C2A_SRR11548778_re-transformed_filtered.csv"
he_complenatery = "/private10/Projects/Gili/HE_workdir/detection/PE_output/filtered_clusters/G2T/G2T_SRR11548778_re-transformed_filtered.csv"

he_data = pd.read_csv(he_filterd_input_file)
he_complementary_data = pd.read_csv(he_complenatery)
he_ids = he_data[["Read_ID", "mate"]]
print(he_ids)
he_complementary_ids = he_complementary_data["Read_ID"].to_list()
fastq_sequences = SeqIO.parse(open(fastq_input_file),'fastq')

flag_loc = 2
count = 0
for record in SeqIO.parse(open(fastq_input_file),'fastq'):
        if record.id in he_ids:
            mate_flag = str(record.description).split()[flag_loc]
            
            if mate_flag == "00":
                  print("both unmapped")
            elif mate_flag == "01":
                  print("mate2 were mapped.")
            else:
                  print("error - mate1 should be umapped")
