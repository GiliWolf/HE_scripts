"""
Author: Gili Wolf
Date: 06-03-2024
Levanon Lab
----------------------------
This script is designed to filter HE read based on pre-defined conditions:
    (-) number of editing sites is bigger than 0
    (-) number of editing sites is bigger than threshold
    (-) Editing fraction bigger than threshold
    (-) phred score of each editing site bigger than threshold
    (-) number of editing sites to read's length ratio
    (-) density of the clusters (length of the cluster to read's length ratio)

the script output 2 CSV file: 
1 - passed csv with the information of the reads that passed all of the consitions 
2 - condition_analysis CSV file with one-hot encoded analysis for each of the conditions
----------------------------

Usage:
    filter_clusters.py -i <SAMPLE_CLUSTERS_CSV> [-o <passed_CSV_PATH> -O <CONDITION_ANALYSIS_PATH>] [Options..]

output files:
    -o --passed_output            path to the read's passed output.
    -O --detailed_condition_output  path to the read's csv file with all the conditions anaylsis output
    -t --output_types"              choose the type of output: 'all', 'passed', 'analysis'

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
from collections import Counter

#Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to filter HE read based on pre-defined conditions.")
# REQUIRED: Input & Output Paths
parser.add_argument("-i", "--input", dest ="sample_clusters_csv", type=str, required=True, help="path to the read's detected cluster (output of detect_clusters.py).")

# OPTIONAL: Conditions parameters
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

# OPTIONAL: control output files
parser.add_argument("-o", "--passed_output", dest ="passed_csv_path", type=str, required=False, default = "passed.csv", help="path to the read's passed output.")
parser.add_argument("-O", "--detailed_condition_output", dest ="condition_analysis_path", type=str, required=False, default = "condition_analysis.csv", help="path to the read's csv file with all the condtions detailed output.")
parser.add_argument("-t", "--output_types", dest="output_types", choices=["all", "passed", "analysis"], default="all", help="choose the type of output: 'all', 'passed', 'analysis'")

args = parser.parse_args()

#TO-DO - 
#   1. 'Edited' column: save reads that were not edited?

# get arguments
sample_clusters_csv = str(args.sample_clusters_csv).strip()
passed_csv_path = str(args.passed_csv_path).strip()
condition_analysis_path = str(args.condition_analysis_path).strip()

# Initialize lists to store rows
passed_reads_data = []
condition_analysis_rows = []

# read data from the detected clusters file
clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)

# output types - 
out_passed = args.output_types == 'all' or args.output_types == 'passed'
out_analysis = args.output_types == 'all' or args.output_types == 'analysis'

# True if consition passed, False if did not
def check_Condition(value, threshold):
    return value >= threshold

def motifs_count(read_seq, editing_sites_map):
    es_positions = list(editing_sites_map.keys())
    nt_counter = {'A': 0, 'G': 0, 'G': 0, 'T':0}
    read_len = len(read_seq)
    counter_dict = {}
    motif_location = ['upstream', 'downstream']
    for loc in motif_location:
        counter_dict[loc] = Counter(nt_counter)
    for pos in es_positions:
        if pos != 0:
            upstream_nt = read_seq[pos - 1]
            counter_dict['upstream'].update([upstream_nt])
        if pos != (read_len - 1):
            downstream_nt = read_seq[pos + 1]
            counter_dict['downstream'].update([downstream_nt])
    data = {}
    for motif_location in counter_dict.keys():
        for nt, count in counter_dict[motif_location].items():
            data[f'{motif_location}_{nt}'] = count
    return data #pd.DataFrame([data])

def es_statistics(editing_sites_map):
    phred_scores = editing_sites_map.values()
    phred_score_avg = sum(phred_scores) / len(phred_scores)
    es_positions = list(editing_sites_map.keys())
    es_positions.sort()
    avg_lambda_func =  lambda es_positions: sum(abs(es_positions[i+1] - es_positions[i]) for i in range(0, len(es_positions) - 1)) / len(es_positions)
    avg_es_distance = avg_lambda_func(es_positions)

    return phred_score_avg, avg_es_distance

    

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
    min_phred_score_passed = not any(edit_site_value < args.min_phred_score for edit_site_value in editing_sites_map.values()) # if one or more is less than the min phred score, not passed
    min_es_length_ratio_passed = check_Condition(read.Number_of_Editing_Sites / read.Alignment_length, args.min_es_length_ratio)
    min_cluster_length_ratio_passed = check_Condition(cluster_len / read.Alignment_length, args.min_cluster_length_ratio)

    #list of conditions -
    conditions_list = [min_editing_sites_passed, min_editing_fraction_passed, min_phred_score_passed, min_es_length_ratio_passed, min_cluster_length_ratio_passed]

    # Update passed_all_conditions flag
    passed_all_conditions = all(conditions_list)

    # if read passed all parameteres add a row to the passed_statisctics
    if passed_all_conditions and out_passed:
        avg_phred_score, avg_es_distance = es_statistics(editing_sites_map)
        motifs_count_data = motifs_count(read.Read_Sequence, editing_sites_map)
        es_statistics_data = {'average_es_phred_score': avg_phred_score, 'average_adjacent_es_distance': avg_es_distance}
        all_data = {}
        all_data.update(es_statistics_data)
        all_data.update(motifs_count_data)
        passed_reads_data.append(all_data)
    # append details to the condition_analysis file
    if out_analysis:
        condition_analysis_rows.append([read[0]] + [passed_all_conditions] + [edited] + conditions_list)

if out_passed:
    # Write all the rows to the passed CSV file
    passed_df = pd.DataFrame(passed_reads_data)
    passed_df.to_csv(passed_csv_path)

if out_analysis:
    # Write all the rows to the condition_analysis CSV file
    with open(condition_analysis_path, 'w', newline='') as condition_analysis_csv:
        csv_writer_condition_analysis = csv.writer(condition_analysis_csv)
        # Write header
        condition_header = ['Read_ID', 'Passed_All', 'Edited', 'Min_Editing_Sites_' + str(args.min_editing_sites), 'Min_Editing_to_Total_MM_Fraction_' + str(args.min_editing_fraction), 'Min_Editing_Phred_Score_' + str(args.min_phred_score), 'Min_Editing_to_Read_Length_Ratio_' + str(args.min_es_length_ratio), 'Min_Cluster_Length_to_Read_Length_Ratio_' + str(args.min_cluster_length_ratio)]
        csv_writer_condition_analysis.writerow(condition_header)
        # Write all rows
        csv_writer_condition_analysis.writerows(condition_analysis_rows)


