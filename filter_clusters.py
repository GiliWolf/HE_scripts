"""
Author: Gili Wolf
Date: 06-03-2024
Levanon Lab
----------------------------
This script is designed to filter HE read based on pre-defined conditions:
    (-) number of editing sites is bigger than 0
    (-) number of editing sites is bigger than threshold
    (-) Editing fracture bigger than threshold
    (-) phred score of each editing site bigger than threshold
    (-) number of editing sites to read's length ratio
    (-) density of the clusters (length of the cluster to read's length ratio)

the script output 2 CSV file: 
1 - filtered csv with the information of the reads that passed all of the consitions 
2 - condition_analysis CSV file with one-hot encoded analysis for each of the conditions
----------------------------

Usage:
    filter_clusters.py -i <SAMPLE_CLUSTERS_CSV> -o <FILTERED_CSV_PATH> -O <CONDITION_ANALYSIS_PATH> [Options..]

optional arguments:
  -es --min_editing_sites           MIN_EDITING_SITES
                                    minimum number of editing sites.

  -ef --min_editing_fraction        MIN_EDITING_FRACTION
                                    minimun fraction of: number of editing sites / total number of MM.

  -ps --min_phred_score             MIN_PHRED_SCORE
                                    minimum value of phred score to each editing site.

  -es2l --min_es_length_ratio       MIN_ES_LENGTH_RATIO
                                    minimum ratio of: number of editing sites to read's length

  -cl2l --min_cluster_length_ratio MIN_CLUSTER_LENGTH_RATIO
                                    minimum ratio of: cluster's length to read's length
"""

import pandas as pd
import csv
import ast
import sys
import argparse
from itertools import combinations

#Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to filter HE read based on pre-defined conditions.")
# Input & Output Paths
parser.add_argument("-i", "--input", dest ="sample_clusters_csv", type=str, required=True, help="path to the read's detected cluster (output of detect_clusters.py).")
parser.add_argument("-o", "--filtered_output", dest ="filtered_csv_path", type=str, required=True, help="path to the read's filtered output.")
parser.add_argument("-O", "--detailed_condition_output", dest ="condition_analysis_path", type=str, required=True, help="path to the read's csv file with all the condtions detailed output.")
# Conditions parameters
# min_editing_sites, default 0 
parser.add_argument("-es", "--min_editing_sites", dest ="min_editing_sites", type=int, required=False, default=0, help="minimum number of editing sites.")
# Editing_to_Total_MM_Fraction, default 0.7 
parser.add_argument("-ef", "--min_editing_fraction", dest ="min_editing_fraction", type=float, required=False, default=0.7, help="minimun fraction of: number of editing sites / total number of MM.")
# min_phred_score, default 30 
parser.add_argument("-ps", "--min_phred_score", dest ="min_phred_score", type=int, required=False, default=30, help="minimum value of phred score to each editing site.")
# min edutung sutes to length ratio, default 0.05 
parser.add_argument("-es2l", "--min_es_length_ratio", dest ="min_es_length_ratio", type=float, required=False, default=0.05, help="minimum ratio of: number of editing sites to read's length")
# ratio of clusters length to read length, default 0.01 
parser.add_argument("-cl2l", "--min_cluster_length_ratio", dest ="min_cluster_length_ratio", type=float, required=False, default=0.1, help="minimum ratio of: cluster's length to read's length")

args = parser.parse_args()

#TO-DO - 
#   2. set of combination parameters  - multiple rung with different parametes. with json?

# get arguments
sample_clusters_csv = str(args.sample_clusters_csv).strip()
filtered_csv_path = str(args.filtered_csv_path).strip()
condition_analysis_path = str(args.condition_analysis_path).strip()

# Initialize lists to store rows
filtered_rows = []
condition_analysis_rows = []

# read data from the detected clusters file
clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)

# True if consition passed, False if did not
def check_Condition(value, threshold):
    return value > threshold

# itearte over each row in the csv file
for read in clusters_df.itertuples():

    # Initialize flag to track if all conditions passed
    passed_all_conditions = True
    edited = False
    # if no editing sites detected - continue to the next read
    # (keep column 'Edited' as False)
    if read.Number_of_Editing_Sites == 0:
        continue
    edited = True
    
    # parse the EditingSites_to_PhredScore_Map (presented as string)
    editing_sites_map = ast.literal_eval(read.EditingSites_to_PhredScore_Map)

    # cluster_len = last position of editing site list minus first position of editing sites list
    cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())

    # Check conditions
    min_editing_sites_passed = check_Condition(read.Number_of_Editing_Sites, args.min_editing_sites)
    min_editing_fraction_passed = check_Condition(read.Editing_to_Total_MM_Fraction, args.min_editing_fraction)
    min_phred_score_passed = any(edit_site_value < args.min_phred_score for edit_site_value in editing_sites_map.values())
    min_es_length_ratio_passed = check_Condition(read.Number_of_Editing_Sites / read.Alignment_length, args.min_es_length_ratio)
    min_cluster_length_ratio_passed = check_Condition(cluster_len / read.Alignment_length, args.min_cluster_length_ratio)

    #list of conditions -
    conditions_list = [min_editing_sites_passed, min_editing_fraction_passed, min_phred_score_passed, min_es_length_ratio_passed, min_cluster_length_ratio_passed]

    # Update passed_all_conditions flag
    passed_all_conditions = all(conditions_list)

    # Append read to filtered rows if all parameteres were passed
    if passed_all_conditions:
        filtered_rows.append(read)
    # append details to the condition_analysis file
    condition_analysis_rows.append([read[0]] + [passed_all_conditions] + [edited] + conditions_list)

# Write all the rows to the filtered CSV file
with open(filtered_csv_path, 'w', newline='') as filtered_csv:
    csv_writer_filtered = csv.writer(filtered_csv)
    # Write header
    csv_writer_filtered.writerow(clusters_df.columns)
    # Write all rows
    csv_writer_filtered.writerows(filtered_rows)

# Write all the rows to the condition_analysis CSV file
with open(condition_analysis_path, 'w', newline='') as condition_analysis_csv:
    csv_writer_condition_analysis = csv.writer(condition_analysis_csv)
    # Write header
    condition_header = ['Read_ID', 'Passed_All', 'Edited' 'Min_Editing_Sites', 'Min_Editing_to_Total_MM_Fraction', 'Min_Editing_Phred_Score', 'Min_Editing_to_Read_Length_Ratio', 'Min_Cluster_Length_to_Read_Length_Ratio']
    csv_writer_condition_analysis.writerow(condition_header)
    # Write all rows
    csv_writer_condition_analysis.writerows(condition_analysis_rows)


