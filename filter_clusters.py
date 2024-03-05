import pandas as pd
import csv
import ast
sample_clusters_csv = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/SRR11548778_A2C_output.csv'

filtered_csv_path = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/filtered_SRR11548778.csv'

rejected_reads_path = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/rejected_SRR11548778.csv'

# Open the Filteres reads output file for writing
filtered_csv = open(filtered_csv_path, "w", newline='')
csv_writer_filtered = csv.writer(filtered_csv)
filtered_header = ['Read_ID', 'Chromosome', 'Position','Alignment_length', 'Reference_Sequence', 'cigar','blocks', 'HE_Sequence', 'Number_of_MM', 'Number_of_Editing_Sites', 'Editing_Fracture', 'Editing_Sites_List']
csv_writer_filtered.writerow(filtered_header)

# Open the Filteres reads output file for writing
rejected_reads = open(rejected_reads_path, "w", newline='')
csv_writer_rejected = csv.writer(rejected_reads)
rejected_header = ['Read_ID', 'Rejected_Condition_List']
csv_writer_rejected.writerow(rejected_header)

clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)
column_names = list(clusters_df.columns)

min_editing_fracture = 0.8
min_phred_score = 30
min_hits_length_ratio = 0.05
min_cluster_length_ratio = 0.1

# itearte over each row in the csv file
for read in clusters_df.itertuples():
    conditions_flag = True
    rejected_condition_list =[]

    if (read.Number_of_Editing_Sites) == 0:
        rejected_condition_list.append('no editing sites')
        csv_writer_rejected.writerow([read[0], rejected_condition_list])
        continue

    editing_sites_map = ast.literal_eval(read.Editing_Sites_List)
    # last position of editing site list minus first position of editing sites list
    cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())
    #FILTER - 
    # ****normalization to the read length??
    # Editing fracture
    if read.Editing_Fracture < min_editing_fracture:
        rejected_condition_list.append('Editing fracture')
        conditions_flag = False

    # phred score of each editing site
    for edit_site_position, edit_site_value in editing_sites_map.items():
        if (edit_site_value < min_phred_score):
            rejected_condition_list.append("phred score")
            conditions_flag = False
            break

    # number of editing sites to read's length ratio
    if (read.Number_of_Editing_Sites / read.Alignment_length) < min_hits_length_ratio:
        rejected_condition_list.append("number of editing sites to read_length ratio")
        conditions_flag = False      

    # density of the clusters:
    if (cluster_len / read.Alignment_length) < min_cluster_length_ratio:
        rejected_condition_list.append("cluster is too dense")
        conditions_flag = False 

    if conditions_flag:
        # Write to output file if all conditions are passed
        csv_writer_filtered.writerow(read)
    else: # write to rejected output file if read did not pass the constions
        csv_writer_rejected.writerow([read[0] , ",".join(condition for condition in rejected_condition_list)])

# Close all the files
filtered_csv.close()
rejected_reads.close()
