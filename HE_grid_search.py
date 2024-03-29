import subprocess
import itertools
import pandas as pd
import os
import json
from concurrent.futures import ThreadPoolExecutor


# Define the parameter grid
param_grid = {
    'min_editing_sites': [0, 1],  
    'min_editing_fraction': [0.5, 0.6], 
    'min_phred_score': [30, 40],  
    'min_es_length_ratio': [0.03, 0.04], 
    'min_cluster_length_ratio': [0.01, 0.02]  
}

# Create combinations of parameters
param_combinations = list(itertools.product(*param_grid.values()))

# outer script
python_command = "python"
filter_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/filter_clusters.py"

# file parameters
input_file = "/private10/Projects/Gili/HE_workdir/detection/second_try/detected_clusters/A2C/A2C_SRR11548778_re-transformed_detected.csv"
sample_id = "A2C_SRR11548778"
output_dir = "/private10/Projects/Gili/HE_workdir/detection/grid_search"


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
max_parallel_threads = 5
with ThreadPoolExecutor(max_workers = max_parallel_threads) as executor:
    executor.map(run_script, param_combinations)


 


