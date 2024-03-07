"""
Author: Gili Wolf
Date: 06-03-2024
Levanon Lab
----------------------------
This script is designed to filter HE read based on pre-defined conditions:
    (-) number of editing sites is bigger than 0
    (-) Editing fracture bigger than threshold
    (-) phred score of each editing site bigger than threshold
    (-) number of editing sites to read's length ratio
    (-) density of the clusters

the script output 2 CSV file: 
first - 
    filtered CSV with the information of the reads which passed all of the conditions,
    with the following parameters:
    1. Read_ID 
    2. Chromosome 
    3. Position
    4. Alignment_length (clipping isn't included, length of the *coding* regions in case of splicing)
    5. Reference_Sequence (clipping isn't included, introns are not incloded, from the genome fasta file)
    6. cigar
    7. blocks (coding regions in case of splicing)
    8. Read_Sequence (original ran read)
    9. Number_of_MM (number of total mismatches of Read_Sequence compare to Reference_Sequence)
    10. Number_of_Editing_Sites (number of mismatches of the ref_base to alt_base kind)
    11. EditingSites_to_PhredScore_Map (key: position, value: phred quality score)

second - 
    rejected CSV with all the reads that didn't pass al least one of the conditions, 
    with a list of the condtions it didn't pass.
----------------------------

Usage:
    python filter_clusters.py <HEreads_csv_path> <output_filtered_csv_path> <output_rejected_reads_path>
"""

import pandas as pd
import csv
import ast
import sys

#TO-DO - 
#   1. parameters can be changed from outside the script - args parse
#   2. set of combination parameters  - multiple rung with different parametes. with json?

# check number of arguments
if len(sys.argv) != 4:
    print("Usage: python filter_clusters.py <HEreads_csv_path> <output_filtered_csv_path> <output_rejected_reads_path>")
    sys.exit(1)

# get arguments
sample_clusters_csv = sys.argv[1]
filtered_csv_path = sys.argv[2]
rejected_reads_path = sys.argv[3]

# Initialize lists to store rows
filtered_rows = []
rejected_rows = []

# Open the Filteres reads output file for writing
filtered_csv = open(filtered_csv_path, "w", newline='')
csv_writer_filtered = csv.writer(filtered_csv)
filtered_header = ['Read_ID', 'Chromosome', 'Position','Alignment_length', 'Reference_Sequence', 'cigar','blocks', 'Read_Sequence', 'Number_of_MM', 'Number_of_Editing_Sites', 'Editing_to_Total_MM_Fraction', 'EditingSites_to_PhredScore_Map']
csv_writer_filtered.writerow(filtered_header)

# Open the Filtered reads output file for writing
rejected_reads = open(rejected_reads_path, "w", newline='')
csv_writer_rejected = csv.writer(rejected_reads)
# Open the rejected reads output file for writing
rejected_header = ['Read_ID', 'Rejected_Condition_List']
csv_writer_rejected.writerow(rejected_header)

# read data from the detected clusters file
clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)

# conditions' parameters
min_editing_sites = 0
min_editing_fracture = 0.8
min_phred_score = 30
min_hits_length_ratio = 0.05
min_cluster_length_ratio = 0.1

# itearte over each row in the csv file
for read in clusters_df.itertuples():
    conditions_flag = True
    rejected_condition_list =[]

    # if no editing sites detectes - continue to the next read
    if (read.Number_of_Editing_Sites) == min_editing_sites:
        rejected_condition_list.append('min editing sites')
        rejected_rows.append([read[0], rejected_condition_list])
        continue
    
    # parse the EditingSites_to_PhredScore_Map (presented as string)
    editing_sites_map = ast.literal_eval(read.EditingSites_to_PhredScore_Map)

    # cluster_len = last position of editing site list minus first position of editing sites list
    cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())
    #FILTER - 
    # ****normalization to the read length??
    # Editing fracture bigger than minimum
    if read.Editing_to_Total_MM_Fraction < min_editing_fracture:
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
        # Append to filtered rows if all conditions are passed
        filtered_rows.append(read)
    else:
        # Append to rejected rows if read did not pass the conditions
        rejected_rows.append([read[0], rejected_condition_list])

# Write all the rows to the filtered CSV file
with open(filtered_csv_path, 'w', newline='') as filtered_csv:
    csv_writer_filtered = csv.writer(filtered_csv)
    # Write header
    csv_writer_filtered.writerow(clusters_df.columns)
    # Write all rows
    csv_writer_filtered.writerows(filtered_rows)

# Write all the rows to the rejected CSV file
with open(rejected_reads_path, 'w', newline='') as rejected_csv:
    csv_writer_rejected = csv.writer(rejected_csv)
    # Write header
    csv_writer_rejected.writerow(['Read_ID', 'Rejected_Condition_List'])
    # Write all rows
    csv_writer_rejected.writerows(rejected_rows)

# Close all the files
filtered_csv.close()
filtered_csv.close()
rejected_reads.close()
