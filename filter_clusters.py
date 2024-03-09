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
import argparse

#Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to filter HE read based on pre-defined conditions.")
# Input & Output Paths
parser.add_argument("-i", "--input", dest ="sample_clusters_csv", type=str, required=True, help="path to the read's detected cluster (output of detect_clusters.py).")
parser.add_argument("-o", "--filtered_output", dest ="filtered_csv_path", type=str, required=True, help="path to the read's filtered output.")
parser.add_argument("-O", "--detailed_condition_output", dest ="rejected_reads_path", type=str, required=True, help="path to the read's csv file with all the condtions detailed output.")
# Conditions parameters
# min_editing_sites, default 0 
parser.add_argument("-es", "--min_editing_sites", dest ="min_editing_sites", type=int, required=False, default=0, help="minimum number of editing sites.")
# Editing_to_Total_MM_Fraction, default 0.7 
parser.add_argument("-ef", "--min_editing_fraction", dest ="min_editing_fraction", type=float, required=False, default=0.7, help="minimun fraction of: number of editing sites / total number of MM.")
# min_phred_score, default 30 
parser.add_argument("-ps", "--min_phred_score", dest ="min_phred_score", type=int, required=False, default=30, help="minimum value of phred score to each editing site.")
# min edutung sutes to length ratio, default 0.05 
parser.add_argument("-es2l", "--min_es_length_ratio", dest ="min_es_length_ratio", type=int, required=False, default=0.05, help="minimum ratio of: number of editing sites to read's length")
# ratio of clusters length to read length, default 0.01 
parser.add_argument("-cl2l", "--min_cluster_length_ratio", dest ="min_cluster_length_ratio", type=float, required=False, default=0.1, help="minimum ratio of: cluster's length to read's length")

args = parser.parse_args()

#TO-DO - 
#   2. set of combination parameters  - multiple rung with different parametes. with json?

# get arguments
sample_clusters_csv = str(args.sample_clusters_csv).strip()
filtered_csv_path = str(args.filtered_csv_path).strip()
rejected_reads_path = str(args.rejected_reads_path).strip()
# conditions' parameters
min_editing_sites = args.min_editing_sites
min_editing_fraction = args.min_editing_fraction
min_phred_score = args.min_phred_score
min_es_length_ratio = args.min_es_length_ratio
min_cluster_length_ratio = args.min_cluster_length_ratio

# Initialize lists to store rows
filtered_rows = []
rejected_rows = []

# read data from the detected clusters file
clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)


def check_Condition(value, condition, threshold ,condition_map):
    # condition passed:
    if value > threshold:
        condition_map[condition] = True
        return condition_map
    # condition didn't pass
    condition_map[condition] = False
    return condition_map


# itearte over each row in the csv file
for read in clusters_df.itertuples():

    # intialize map to store whether a condition is passed or not :
    condition_map = {'Edited': False, 'Min_Editing_Sites': False, 'Min_Editing_to_Total_MM_Fraction': False, 'Min_Editing_Phred_Score': False, 'Min_Editing_to_Read_Length_Ratio': False, 'Min_Cluster_Length_to_Read_Length_Ratio': False, 'Passed_All': False}

    # if no editing sites detectes - continue to the next read
    if (read.Number_of_Editing_Sites) == 0:
        continue
    else:
        condition_map['Edited'] = True
    
    # parse the EditingSites_to_PhredScore_Map (presented as string)
    editing_sites_map = ast.literal_eval(read.EditingSites_to_PhredScore_Map)

    # cluster_len = last position of editing site list minus first position of editing sites list
    cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())


    #FILTER - 
    # Number of editing sites bigger than minimum
    condition_map = check_Condition(read.Number_of_Editing_Sites,
                                    'Min_Editing_Sites',
                                    min_editing_sites, 
                                    condition_map)
    
    # Editing fracture bigger than minimum
    condition_map = check_Condition(read.Editing_to_Total_MM_Fraction,
                                    'Min_Editing_to_Total_MM_Fraction',
                                    min_editing_fraction, 
                                    condition_map)
    
    # phred score of each editing sites
    for edit_site_position, edit_site_value in editing_sites_map.items():
        condition_map = check_Condition(edit_site_value,
                                        'Min_Editing_Phred_Score',
                                        min_phred_score,
                                        condition_map)
    
    # number of editing sites to read's length ratio
    condition_map = check_Condition(read.Number_of_Editing_Sites / read.Alignment_length,
                                    'Min_Editing_to_Read_Length_Ratio',
                                    min_es_length_ratio,
                                    condition_map)

    # density of the clusters:
    condition_map = check_Condition(cluster_len / read.Alignment_length,
                                    'Min_Cluster_Length_to_Read_Length_Ratio',
                                    min_cluster_length_ratio,
                                    condition_map)

    if all(condition_map.values()):
        # Append to filtered rows if all conditions are passed
        condition_map['Passed_All'] = True
        filtered_rows.append(read)
    # [passed_flag for passed_flag in condition_map.values()]
    rejected_rows.append([passed_flag for passed_flag in condition_map.values()])


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
    csv_writer_rejected.writerow(list(condition_map.keys()))
    # Write all rows
    csv_writer_rejected.writerows(rejected_rows)


