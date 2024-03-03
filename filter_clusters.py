import pandas as pd

sample_clusters_csv = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/SRR11548778_A2C_output.csv'

filtered_csv_path = '/private10/Projects/Gili/HE_workdir/detection/detect_clusters_test/filtered_SRR11548778.csv'

clusters_df = pd.read_csv(sample_clusters_csv, index_col=0)
column_names = list(clusters_df.columns)

min_editing_fracture = 0.6
min_phred_score = 30
min_hits_length_ratio = 0.05
min_cluster_length_ratio = 0.1

# itearte over each row in the csv file
for read in clusters_df.itertuples():
    editing_sites_map = map(read.Editing_Sites_List)
    cluster_len = max(editing_sites_map.keys()) - min(editing_sites_map.keys())
    #FILTER - 
    # normalization to the read length??
    # Editing fracture
    assert(read.Editing_Fracture >= min_editing_fracture)

    # phred score of each editing site
    for edit_site in map(read.Editing_Sites_List):
        assert(edit_site.value >= min_phred_score)

    # number of editing sites to read's length ratio
    assert((read.Number_of_Editing_Sites / read.Alignment_length) >= min_hits_length_ratio)

    # density of the clusters:
    assert((cluster_len / read.Alignment_length) >= min_cluster_length_ratio)