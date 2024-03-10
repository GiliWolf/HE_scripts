import subprocess
import itertools
import pandas as pd
import os
import json

# for creating json files:
# import json

# data = {'name': 'John', 'age': 30}

# with open('data.json', 'w') as f:
#     json.dump(data, f)


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

python_command = "python"
filter_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/filter_clusters.py"

input_file = "/private10/Projects/Gili/HE_workdir/detection/second_try/detected_clusters/A2C/A2C_SRR11548778_re-transformed_detected.csv"
sample_id = "A2C_SRR11548778"
output_dir = "/private10/Projects/Gili/HE_workdir/detection/grid_search"

#  Run the script for each parameter combination
results = []
for params in param_combinations:
    # Construct the command to run the script with the current parameters
    command = [python_command, filter_script]
    command_parameters = []
    path_parameters=''
    for param, value in zip(param_grid.keys(), params):
        str_value = str(value)
        path_parameters+= '_'
        path_parameters+= str_value
        command_parameters.extend(['--' + param, str_value])
    filtered_output_path = os.path.join(output_dir, sample_id + path_parameters + "_filtered")
    condition_analysis_output_path = os.path.join(output_dir, sample_id + path_parameters + "_condition_analysis")

    command.extend(['-i', input_file])
    command.extend(['-o', filtered_output_path])
    command.extend(['-O', condition_analysis_output_path])
    command.extend(command_parameters)

    def analyze_Files():
        df = pd.read_csv('/private10/Projects/Gili/HE_workdir/detection/grid_search/A2C_SRR11548778_0_0.5_30_0.03_0.01_condition_analysis')
        # Get the list of condition columns
        condition_columns = df.columns[2:]

        json_data = {}
        number_of_filtered_reads = df['Passed_All'].sum()
        json_data["number of filtererd reads: "] = int(number_of_filtered_reads)

        each_condition = df.drop(columns=['Read_ID', 'Passed_All']).sum()
        rc_condition_map ={}
        for condition, reads_count in each_condition.items():
            rc_condition_map[condition] = int(reads_count)
        json_data["number of filtererd reads for each condition:"] = rc_condition_map

        # Create a dictionary to store the number of reads that passed all conditions for each subset
        subset_passed_counts = {}
        m = len(condition_columns) - 1
        # Iterate over all possible subsets of M conditions
        for subset in itertools.combinations(condition_columns, m):
            # Calculate the number of reads that passed all conditions in the current subset
            subset_key = "_".join(subset)
            passed_subset = df[df[list(subset)].all(axis=1)]
            num_passed_subset = passed_subset.shape[0]
            # Store the count in the dictionary
            subset_passed_counts[subset_key] = int(num_passed_subset)
        json_data["Number of reads that passed each subset of conditions:"] = subset_passed_counts

        with open('data.json', 'w') as f:
            json.dump(json_data, f, indent = 4)

    analyze_Files()
    break

    # try:
    #     process = subprocess.run(command, capture_output=True, text=True, check=True)
    #     analyze_Files()
    # # Check if the process completed successfully
    #     if process.returncode != 0:
    #         print("Process failed with return code:", process.returncode)
    #     # # Print stdout
    #     # print("STDOUT:")
    #     # print(process.stdout)
    #     # # Print stderr
    #     # print("STDERR:")
    #     # print(process.stderr)
    # except subprocess.CalledProcessError as e:
    #     # Print error message if the process failed
    #     print("Error:", e)
    #     print("STDERR from subprocess:")
    #     print(e.stderr)

    # break