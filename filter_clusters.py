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

the script output 3 files (if -t == 'all'): 
1 - CSV with the information of the reads that passed all of the consitions 
2 - CSV file with True/False for each of the conditions (passes/not passed) for each read
3 - JSON summary file of the sample
----------------------------

Usage:
    filter_clusters.py -i <SAMPLE_CLUSTERS_CSV> \
                       [-o <passed_CSV_PATH> -O <CONDITION_ANALYSIS_PATH> -id <SAMPLE_ID> \
                        -j <JSON_PATH> -t {all,passed,analysis,summary}] \
                       [Parameters Options..]

Arguments: 
REQUIRED-
input files:
    -i --input                      path to the read's detected clusters (output of detect_clusters.py)

OPTIONAL-
output files:
    -id --sample_id                 sample ID
    -o --passed_output              path to the read's passed output.
    -O --detailed_condition_output  path to the read's csv file with all the conditions anaylsis output
    -j --json_summary               path to the sample's summary json file
    -t --output_types               choose the type of output: 'all', 'passed', 'analysis', 'summary'
filter parameters:
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
import json


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
parser.add_argument("-id", "--sample_id",dest = "sample_id", required=False, type=str, default="sample", help="Sample ID.")
parser.add_argument("-o", "--passed_output", dest ="passed_csv_path", type=str, required=False, default = "passed.csv", help="path to the read's passed output.")
parser.add_argument("-O", "--detailed_condition_output", dest ="condition_analysis_path", type=str, required=False, default = "condition_analysis.csv", help="path to the read's csv file with all the condtions detailed output.")
parser.add_argument("-j", "--json_summary", dest="json_path", type=str, required=False,default = "summary.json", help="path to the sample's summary json file.")
parser.add_argument("-t", "--output_types", dest="output_types", choices=["all", "passed", "analysis", "summary"], default="all", help="choose the type of output: 'all', 'passed', 'analysis', 'summary'")


args = parser.parse_args()

#TO-DO - 
#   1. 'Edited' column: save reads that were not edited?

# get arguments
sample_clusters_csv = str(args.sample_clusters_csv).strip()
passed_csv_path = str(args.passed_csv_path).strip()
condition_analysis_path = str(args.condition_analysis_path).strip()
json_summary_path = str(args.json_path).strip()

# Initialize lists to store rows
passed_reads_data = []
condition_analysis_rows = []

# read data from the detected clusters file
clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)

# output types -
sample_id = args.sample_id
out_passed = args.output_types == 'all' or args.output_types == 'passed'
out_analysis = args.output_types == 'all' or args.output_types == 'analysis'
out_json = args.output_types == 'all' or args.output_types == 'summary'
out_json = True
# script parameters 
motif_location = ['upstream', 'downstream']

# True if consition passed, False if did not
def check_Condition(value, threshold):
    return value >= threshold

def sample_motif(nt_count_df, location, total_num_of_ES):
    nt_precentage = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    for nt in nt_precentage.keys():
        perc = (nt_count_df[f'{location}_{nt}'].sum() / (total_num_of_ES))
        nt_precentage[nt] = perc
    return nt_precentage

#create json file summarizing key statistics of he passed rfeads
def create_json(passed_df, condition_df):

    # init general data dictonary 
    json_data = {"sample_id: ": sample_id}
    
    # add total number of filtered reads
    json_data["number of passed reads: "] = int(condition_df['Passed_All'].sum())

    # average es number
    average_num_of_ES = passed_df['number_of_ES'].mean()
    json_data["average number of ES: "] = average_num_of_ES

    # get sample ES motif
    total_num_of_ES = passed_df['number_of_ES'].sum()
    for loc in motif_location:
        nt_count_df = passed_df[[f'{loc}_{nt}' for nt in ['A', 'C', 'G', 'T']]]
        nt_precentage =sample_motif(nt_count_df, loc, total_num_of_ES)
        json_data[f'{loc}_motif: '] = nt_precentage

    # Get the list of condition columns
    condition_df = condition_df.drop(columns=['Read_ID', 'Passed_All'])

    # get values of parameters for each column and remove it from the name
    parameters_values = {}
    for col in condition_df.columns:
        if (col == 'Edited'):
            continue
        # extract param's value
        param_value = col.rsplit('_', 1)[-1]
        # rename the columns to discard the param's value
        column_name = col.rsplit('_', 1)[0]
        condition_df.rename(columns={col: column_name}, inplace=True)
        # add to the parameter's values
        parameters_values[column_name] = param_value
    json_data["Parameters' Values: "] = parameters_values

    # Create dictionary with the sum of reads which passed each condition
    each_condition = condition_df.sum().astype(int).to_dict()
    json_data["number of passed reads for each condition"] = each_condition

    return json_data


# count number of upstream + downstrean nt for each ES, returns a dictionary with the sum for each type of motif nt 
def motifs_count(read_seq, editing_sites_map):
    es_positions = list(editing_sites_map.keys())
    nt_counter = {'A': 0, 'C': 0, 'G': 0, 'T':0}
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
    return data 

# calculat and returns average phred score for all of the ES, and averge of al distances netween adjactent ES
def es_statistics(editing_sites_map):
    phred_scores = editing_sites_map.values()
    phred_score_avg = sum(phred_scores) / len(phred_scores)
    es_positions = list(editing_sites_map.keys())
    es_positions.sort()
    avg_lambda_func =  lambda es_positions: sum(abs(es_positions[i+1] - es_positions[i]) for i in range(0, len(es_positions) - 1)) / len(es_positions)
    avg_es_distance = avg_lambda_func(es_positions)
    num_of_es = len(editing_sites_map)
    return phred_score_avg, avg_es_distance, num_of_es

    

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
        avg_phred_score, avg_es_distance, num_of_es = es_statistics(editing_sites_map)
        motifs_count_data = motifs_count(read.Read_Sequence, editing_sites_map)
        es_statistics_data = {'number_of_ES': num_of_es, 'average_es_phred_score': avg_phred_score, 'average_adjacent_es_distance': avg_es_distance}
        all_data = {"Read_ID": read[0]}
        all_data.update(es_statistics_data)
        all_data.update(motifs_count_data)
        passed_reads_data.append(all_data)
    # append details to the condition_analysis file
    if out_analysis:
        condition_analysis_rows.append([read[0]] + [passed_all_conditions] + [edited] + conditions_list)


# WRITE OUTPUT
#CREATE DATA FRAMES
passed_df = pd.DataFrame(passed_reads_data)

condition_header = ['Read_ID', 'Passed_All', 'Edited', 'Min_Editing_Sites_' + str(args.min_editing_sites), 'Min_Editing_to_Total_MM_Fraction_' + str(args.min_editing_fraction), 'Min_Editing_Phred_Score_' + str(args.min_phred_score), 'Min_Editing_to_Read_Length_Ratio_' + str(args.min_es_length_ratio), 'Min_Cluster_Length_to_Read_Length_Ratio_' + str(args.min_cluster_length_ratio)]
condition_df = pd.DataFrame(condition_analysis_rows, columns=condition_header)

if out_passed:
    # Write all the rows to the passed CSV file
    passed_df.to_csv(passed_csv_path)
if out_analysis:
    # Write all the rows to the condition_analysis CSV file
    condition_df.to_csv(condition_analysis_path)
if out_json:
    json_data = create_json(passed_df, condition_df)
    # write the data to a json file
    with open(json_summary_path, 'w') as f:
        json.dump(json_data, f, indent = 4)


