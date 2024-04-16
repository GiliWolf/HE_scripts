"""
Author: Gili Wolf
Date: 30-03-2024
Levanon Lab
----------------------------
Description:
    This script is designed to parallel run the HE_filter_clusters script for all combinations of given parameter ranges.

Algorithm:
    1. Extract combinations of parameters from the grid.
    2. Run the script using these parameters.
    3. Create a statistics JSON file.
    4. Delete raw output files.

For each set of parameters, the script outputs a JSON file containing:
    - Number of filtered reads
    - Parameters' values of the run
    - Number of filtered reads for each condition
----------------------------

usage: 
    HE_grid_search.py [-h] -i INPUT_FILE -id SAMPLE_ID -O OUTPUT_DIR [-es X X X] [-ef X X X] [-ps X X X] [-es2l X X X] [-cl2l X X X]

REQUIRED ARGUMENTS:
  -i INPUT_FILE, --input INPUT_FILE
                        Path to the input CSV file.
  -id SAMPLE_ID, --sample_id SAMPLE_ID
                        Sample ID.
  -D OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to the output directory.
OPTIONAL ARGUMENTS:
  -th MAX_THREADS, --threads MAX_THREADS
                        max number of threads. (default - 5)
  -es X X X, --min_editing_sites X X X
                        range for minimum number of editing sites.
  -ef X X X, --min_editing_fraction X X X
                        range for minimum fraction of: number of editing sites / total number of MM.
  -ps X X X, --min_phred_score X X X
                        range for minimum value of phred score to each editing site.
  -es2l X X X, --min_es_length_ratio X X X
                        range for minimum ratio of: number of editing sites to read's length
  -cl2l X X X, --min_cluster_length_ratio X X X
                        range for minimum ratio of: cluster's length to read's length
"""
import subprocess
import itertools
import pandas as pd
import os
import json
import sys
from concurrent.futures import ThreadPoolExecutor
import argparse
import numpy as np
# TO-DO:
#   2- decieded which files to output
#   3- add summarizing file

# Controling Arguments:
parser = argparse.ArgumentParser(description="""
                                 This script is designed to grid-search HE filtered reads statitics, for each combination of parameters.
                                 for each arguments you need to provide min_value, max_value and size of interval.
                                 """)
# IO FILES
parser.add_argument("-i","--input",dest = "input_file", required=True, type=str, help="Path to the input CSV file.")
parser.add_argument("-id", "--sample_id",dest = "sample_id", required=True, type=str, help="Sample ID.")
parser.add_argument("-D", "--output_dir", dest = "output_dir", required=True, type=str, help="Path to the output directory.")
# PARALLEL PARAMETERS:
parser.add_argument("-th", "--threads", dest = "max_threads", required=False, type=int, default=5, help="max number of threads.")

# Define Grid arguments with min, max, and interval size
arg_infos = [
    ("-es", "--min_editing_sites", "min_editing_sites", int, 0, 1, 1, "range for minimum number of editing sites."),
    ("-ef", "--min_editing_fraction", "min_editing_fraction", float, 0.5, 0.6, 0.1, "range for minimum fraction of: number of editing sites / total number of MM."),
    ("-ps", "--min_phred_score", "min_phred_score", int, 30, 40, 10, "range for minimum value of phred score to each editing site."),
    ("-es2l", "--min_es_length_ratio", "min_es_length_ratio", float, 0.03, 0.04, 0.01, "range for minimum ratio of: number of editing sites to read's length"),
    ("-cl2l", "--min_cluster_length_ratio", "min_cluster_length_ratio", float, 0.01, 0.02, 0.01, "range for minimum ratio of: cluster's length to read's length")
]
metavar_symbol = 'X'
for arg_short, arg_long, dest, arg_type, min_val, max_val, step, help_msg in arg_infos:
    parser.add_argument(arg_short, arg_long, dest=dest, type=arg_type, nargs=3, required=False, default=[min_val, max_val, step],metavar=metavar_symbol, help=help_msg)

args = parser.parse_args()

# Generate parameter grid
param_grid = {}
for arg_short, arg_long, dest, arg_type, min_val, max_val, step, help in arg_infos:
    min_val, max_val, step = getattr(args, dest)
    # using step/2 upper limit in order to prevent float unexpected behavor (from numpy docs: "where step is not an integer and floating point round-off affects the length of out")
    param_grid[dest] = list(np.arange(min_val, max_val + (step/2), step, dtype=arg_type))

param_combinations = list(itertools.product(*param_grid.values()))

# outer script
python_command = "python"
filter_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/filter_clusters.py"

# file parameters
# input_file = "/private10/Projects/Gili/HE_workdir/detection/second_try/detected_clusters/A2C/A2C_SRR11548778_re-transformed_detected.csv"
# sample_id = "A2C_SRR11548778"
# output_dir = "/private10/Projects/Gili/HE_workdir/detection/grid_search"
input_file = args.input_file
sample_id = args.sample_id
output_dir = args.output_dir
output_files_type = 'analysis'


#create json file summarizing key statistics from the output file
def create_json(analysis_file, json_path):
    df = pd.read_csv(analysis_file)

    # init general data dictonary 
    json_data = {}

    # add total number of filtered reads
    json_data["number of filtererd reads: "] = int(df['Passed_All'].sum())

    # Get the list of condition columns
    conditions_df = df.drop(columns=['Read_ID', 'Passed_All'])

    # get values of parameters for each column and remove it from the name
    parameters_values = {}
    for col in conditions_df.columns:
        if (col == 'Edited'):
            continue
        # extract param's value
        param_value = col.rsplit('_', 1)[-1]
        # rename the columns to discard the param's value
        column_name = col.rsplit('_', 1)[0]
        conditions_df.rename(columns={col: column_name}, inplace=True)
        # add to the parameter's values
        parameters_values[column_name] = param_value
    json_data["Parameters' Values: "] = parameters_values

    # Create dictionary with the sum of reads which passed each condition
    each_condition = conditions_df.sum().astype(int).to_dict()
    json_data["number of filtered reads for each condition"] = each_condition

    # # Create a dictionary to store the number of reads that passed all conditions for each subset
    # subset_passed_counts = {}
    # m = len(conditions_df) - 1
    # # Iterate over all possible subsets of M conditions
    # for subset in itertools.combinations(conditions_df, m):
    #     # Calculate the number of reads that passed all conditions in the current subset
    #     subset_key = "_".join(subset)
    #     passed_subset = df[df[list(subset)].all(axis=1)]
    #     num_passed_subset = passed_subset.shape[0]
    #     # Store the count in the dictionary
    #     subset_passed_counts[subset_key] = int(num_passed_subset)
    # json_data["Number of reads that passed each subset of conditions:"] = subset_passed_counts

    # write the data to a json file
    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent = 4)


def run_script(params):
    # Construct the command to run the script with the current parameters
    command = [python_command, filter_script]
    command_parameters = []
    path_parameters=''
    for param, value in zip(param_grid.keys(), params):
        str_value = str(value)
        # add paramter value to the output file name
        path_parameters+= ('_' + str_value)
        # add parameter to the command
        command_parameters.extend(['--' + param, str_value])
    # construct output path files
    filtered_output_path = os.path.join(output_dir, sample_id + "_filtered" + path_parameters + ".csv")
    condition_analysis_output_path = os.path.join(output_dir, sample_id + "_condition_analysis" + path_parameters + ".csv")
    json_output_path = os.path.join(output_dir, sample_id + path_parameters + ".json")
    # concat all parameters: input, output, grid search
    command.extend(['-i', input_file])
    command.extend(['-o', filtered_output_path])
    command.extend(['-O', condition_analysis_output_path])
    command.extend(['-t', output_files_type])
    command.extend(command_parameters)

    try:
        # run process
        process = subprocess.run(command, capture_output=True, text=True, check=True)
        # create json
        create_json(condition_analysis_output_path,json_output_path)
        # remove file after json analysis 
        os.remove(condition_analysis_output_path)
    # Check if the process completed successfully
        if process.returncode != 0:
            print("Process failed with return code:", process.returncode)
        # Print stdout
        print("STDOUT:")
        print(process.stdout)
        # Print stderr
        print("STDERR:")
        print(process.stderr)
    except subprocess.CalledProcessError as e:
        # Print error message if the process failed
        print("Error:", e)
        print("STDERR from subprocess:")
        print(e.stderr)
    except FileNotFoundError:
        print(f"File {condition_analysis_output_path} does not exist.")


# Run the script for each parameter combination concurrently
with ThreadPoolExecutor(max_workers = args.max_threads) as executor:
    executor.map(run_script, param_combinations)


 


