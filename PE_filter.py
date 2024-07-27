from Bio import SeqIO
import pandas as pd

"""
to-do:
      1) get info from detections file 
      2) check resonable insert size 
      3) check reversed editing ** only check this, add column 
      4) if 1 is mapped - check position are ok
"""


fastq_input_file = "/private10/Projects/Gili/HE_workdir/first_part/PE_test/first_map/SRR11548778.Unmapped.out.mate1"
he_filterd_input_file = "/private10/Projects/Gili/HE_workdir/detection/PE_output/filtered_clusters/C2A/C2A_SRR11548778_re-transformed_filtered.csv"
he_complenatery = "/private10/Projects/Gili/HE_workdir/detection/PE_output/filtered_clusters/G2T/G2T_SRR11548778_re-transformed_filtered.csv"

he_data = pd.read_csv(he_filterd_input_file)
he_complementary_data = pd.read_csv(he_complenatery)
he_ids = he_data[["Read_ID", "mate"]]
print(he_ids)
he_complementary_ids = he_complementary_data["Read_ID"].to_list()
fastq_records = SeqIO.parse(open(fastq_input_file),'fastq')
fastq_dict = {record.id: record for record in fastq_records}
flag_loc = 2
count = 0


def both_unmapped(read_id):
      mate1 = he_data[he_data["Read_ID"] == read_id]
      mate2 = he_complementary_data[he_complementary_data["Read_ID"] == read_id]
      print(f"mate1: {mate1}\nmate2: {mate2}")


for row in he_ids.itertuples():
      read_id = row.Read_ID
      record = fastq_dict[read_id]
      mate_flag = str(record.description).split()[flag_loc]
            
      if mate_flag == "00":
            if read_id in he_complementary_ids:
                  both_unmapped(read_id)
      elif mate_flag == "01":
            print("mate2 were mapped.")
      else:
            print("error - mate1 should be umapped")

    

