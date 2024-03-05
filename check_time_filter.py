import timeit

# Define setup code
setup_code = """
import pandas as pd

sample_clusters_csv = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/SRR11548778_A2C_output.csv'
filtered_csv_path = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/filtered_SRR11548778.csv'

clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)
column_names = list(clusters_df.columns)

min_editing_fracture = 0.6
min_phred_score = 30
min_hits_length_ratio = 0.05
min_cluster_length_ratio = 0.1
"""

# Define code to be timed
code_if_approach = """
with open(filtered_csv_path, 'w') as output_file:
    for read in clusters_df.itertuples():
        editing_sites_map = map(read.Editing_Sites_List)
        cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())
        conditions_met = True
        if not (read.Editing_Fracture >= min_editing_fracture):
            conditions_met = False
        for edit_site in map(read.Editing_Sites_List):
            if not (edit_site.value >= min_phred_score):
                conditions_met = False
        if not ((read.Number_of_Editing_Sites / read.Alignment_length) >= min_hits_length_ratio):
            conditions_met = False
        if not ((cluster_len / read.Alignment_length) >= min_cluster_length_ratio):
            conditions_met = False
        if conditions_met:
            output_file.write(','.join(str(field) for field in read))
"""

code_all_approach = """
with open(filtered_csv_path, 'w') as output_file:
    for read in clusters_df.itertuples():
        editing_sites_map = map(read.Editing_Sites_List)
        cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())
        conditions_met = all([
            read.Editing_Fracture >= min_editing_fracture,
            all(edit_site.value >= min_phred_score for edit_site in map(read.Editing_Sites_List)),
            (read.Number_of_Editing_Sites / read.Alignment_length) >= min_hits_length_ratio,
            (cluster_len / read.Alignment_length) >= min_cluster_length_ratio
        ])
        if conditions_met:
            output_file.write(','.join(str(field) for field in read))
"""

# Measure execution time
time_if_approach = timeit.timeit(stmt=code_if_approach, setup=setup_code, number=100)
time_all_approach = timeit.timeit(stmt=code_all_approach, setup=setup_code, number=100)

# Print results
print("Execution time for if approach:", time_if_approach)
print("Execution time for all approach:", time_all_approach)
