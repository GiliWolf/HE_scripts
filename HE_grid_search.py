import subprocess
import itertools
import pandas as pd
import os

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
        df = pd.read_csv(filtered_output_path)
        number_of_filtered_reads = df.shape[0]
        print("number of filtererd reads: ", number_of_filtered_reads)

    # run
    # process = subprocess.run(command, capture_output=True, text=True)
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=True)
        analyze_Files()
    # Check if the process completed successfully
        if process.returncode != 0:
            print("Process failed with return code:", process.returncode)
        # # Print stdout
        # print("STDOUT:")
        # print(process.stdout)
        # # Print stderr
        # print("STDERR:")
        # print(process.stderr)
    except subprocess.CalledProcessError as e:
        # Print error message if the process failed
        print("Error:", e)
        print("STDERR from subprocess:")
        print(e.stderr)

    break