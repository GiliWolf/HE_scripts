from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import ast
import sys
import argparse
from collections import Counter
import json
import numpy as np


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

the script output 3 files (if -t == 'all'): 
1 - CSV with the information of the reads that passed all of the consitions 
2 - CSV file with True/False for each of the conditions (passes/not passed) for each read
3 - JSON summary file of the sample
----------------------------
Multimappers:
Before applying the filtration process, if multiple loci exist for a given sample ID, the read with the maximum editing fraction will be selected. The editing fraction is calculated based on all MM, without considering quality checks.
----------------------------
Parallel Processing:
The CSV file will be divided into batches for parallel processing based on the number of threads. The size of each batch is determined as follows: int(num_of_reads / max_threads).
If -b is set, will use the number specified as batch size.
----------------------------

Usage:
    filter_clusters.py -i <SAMPLE_CLUSTERS_CSV> \
                       [-o <passed_CSV_PATH> -O <CONDITION_ANALYSIS_PATH> -id <SAMPLE_ID> \
                        -j <JSON_PATH> -t {all,passed,analysis,summary} -t <MAX_THREADS> -b <BATCH_SIZE>] \
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
    -ot --output_types              choose the type of output: 'all', 'passed', 'analysis', 'summary'
    -t --threads                    Number of threads to parallel processing.(default: 5)
    -b --batch                      Batch size (number of reads) to be processed by each thread. (default: 50)
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
#Controling Arguments:
parser = argparse.ArgumentParser(description="This script is designed to filter HE read based on pre-defined conditions.")
# REQUIRED: Input & Output Paths
parser.add_argument("-i", "--input", dest ="sample_clusters_csv", type=str, required=True, help="path to the read's detected cluster (output of detect_clusters.py).")

# OPTIONAL: Controlling parallel processing
parser.add_argument("-t", "--threads", type=int, required=False, default=5, help="Number of threads to parallel processing.(default: 5)")
parser.add_argument("-b", "--batch", type=int, required=False, default=50, help="Specify the batch size, i.e., the number of reads to be processed by each thread. If not set, will use: int(num_of_reads / max_threads) (default: 50)")


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
parser.add_argument("-M", "--motifs_output", dest ="motifs_analysis_path", type=str, required=False, default = "motifs_analysis.csv", help="path to the output csv file with each read's motifs (count of each nucleotide before and after ES).")
parser.add_argument("-B", "--bed_output", dest ="BED_clusters_path", type=str, required=False, default = "HE_clusters.bed", help="path to the output BED file.")
parser.add_argument("-j", "--json_summary", dest="json_path", type=str, required=False,default = "summary.json", help="path to the sample's summary json file.")
parser.add_argument("-ot", "--output_types", dest="output_types", choices=["all", "passed", "analysis","motifs", "bed", "summary"], default="all", help="choose the type of output: 'all', 'passed', 'analysis', 'motifs', 'bed',  'summary'")

# get arguments
args = parser.parse_args()
sample_clusters_csv = str(args.sample_clusters_csv).strip()
passed_csv_path = str(args.passed_csv_path).strip()
condition_analysis_path = str(args.condition_analysis_path).strip()
motifs_csv_path =str(args.motifs_analysis_path).strip()
bed_path=str(args.BED_clusters_path).strip()
json_summary_path = str(args.json_path).strip()

# parallel processing controling
max_threads = args.threads
user_batch_size = args.batch

# output types -
sample_id = args.sample_id
out_passed = args.output_types == 'all' or args.output_types == 'passed'
out_analysis = args.output_types == 'all' or args.output_types == 'analysis'
out_json = args.output_types == 'all' or args.output_types == 'summary'
out_motifs = args.output_types == 'all' or args.output_types == 'motifs'
out_bed = args.output_types == 'all' or args.output_types == 'bed'

# INIT SCRIPT PARAMETERS
# Initialize lists to store rows
passed_reads_data = []
motifs_data= []
bed_file_data =[]
condition_analysis_rows = []
motif_location = ['upstream', 'downstream']


# HELPER FUNCIONS
# True if consition passed, False if did not
def check_Condition(value, threshold):
    return value >= threshold

def sample_motif(nt_count_df, location, total_num_of_ES):
    nt_precentage = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    nt_count = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    for nt in nt_precentage.keys():
        total_nt_count = int(nt_count_df[f'{location}_{nt}'].sum())
        perc = (total_nt_count / total_num_of_ES)
        nt_precentage[nt] = perc
        nt_count[nt] = total_nt_count
    return nt_precentage, nt_count

# create json file summarizing key statistics of HE passed reads
def create_json(passed_df, motifs_df, condition_df,multimappers_count, unique_count):

    # init general data dictonary 
    json_data = {"sample_id:": sample_id,"multimappers total reads before filtering:": multimappers_count, "uniquely total reads before filtering:": unique_count}
    
    # add total number of filtered reads
    json_data["number of passed reads:"] = int(condition_df['Passed_All'].sum())

    # average es number
    average_num_of_ES = passed_df['Number_of_Passed_ES'].mean()
    json_data["average number of ES:"] = average_num_of_ES

    # get sample ES motif
    total_num_of_ES = passed_df['Number_of_Passed_ES'].sum()
    for loc in motif_location:
        nt_count_df = motifs_df[[f'{loc}_{nt}' for nt in ['A', 'C', 'G', 'T']]]
        nt_precentage, nt_count = sample_motif(nt_count_df, loc, total_num_of_ES)
        json_data[f'{loc}_motif_count:'] = nt_count
        json_data[f'{loc}_motif_percentage:'] = nt_precentage

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
    json_data["Parameters' Values:"] = parameters_values

    # Create dictionary with the sum of reads which *did not* pass each condition
    json_data["Filtering Results (did not pass the condition):"] = (condition_df == False).sum().astype(int).to_dict()

    return json_data

# create a bed line of the HE clusters (the genomic interval of the closest editing sites)
def parse_bed_row(read, passed_es_map):
    genomic_blocks = ast.literal_eval(read.Genomic_Position_Splicing_Blocks_0based)
    read_blocks = ast.literal_eval(read.Read_Relative_Splicing_Blocks_0based)
    rows = []
    for i, block in enumerate(read_blocks):
        block_start, block_end = block[0], block[1]
        genomic_start = genomic_blocks[i][0]

        # find cluster: first es in the block, last es in the block, and number of es in the block
        first_fit_pos = None  # Initialize to None to check later
        last_fit_pos = None
        cluster_score = 0
        for pos in passed_es_map:
            if block_start <= pos <= block_end:
                cluster_score += 1  # Increment the sum for each matching position
                if first_fit_pos is None:  # If first_fit_pos hasn't been set yet
                    first_fit_pos = pos
                last_fit_pos = pos  # update last_fit_pos to the current match
        print(f"block ${i}: first_fit_pos = ${first_fit_pos}")
        print(f"block ${i}: last_fit_pos = ${last_fit_pos}")

        if cluster_score==0:
            continue
        rows.append({
            'chr': read.Chromosome,
            'start': genomic_start + first_fit_pos - block_start,
            'end': genomic_start + last_fit_pos - block_start,
            'name': read.Read_ID,
            'score': cluster_score,
            'strand':read.Strand
        })  
    return rows
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

# calculate and returns average phred score for all of the ES, and averge of al distances between adjactent ES
def es_statistics(editing_sites_map):
    phred_scores = editing_sites_map.values()
    phred_score_avg = sum(phred_scores) / len(phred_scores)
    es_positions = list(editing_sites_map.keys())
    es_positions.sort()
    avg_lambda_func =  lambda es_positions: sum(abs(es_positions[i+1] - es_positions[i]) for i in range(0, len(es_positions) - 1)) / len(es_positions)
    avg_es_distance = avg_lambda_func(es_positions)
    return phred_score_avg, avg_es_distance

# for each site checks if the phred score is bigger then the threshold
# if so  - consider it as passed site 
# returns number of passed sites and the position-phred score map of the passed sites
def filter_sites(sites_map, score_threshold):
    passed_sites_map = {pos: score for pos, score in sites_map.items() if score >= score_threshold}
    return len(passed_sites_map), passed_sites_map

def process_reads(clusters_df):
    # itearte over each row in the csv file
    for read in clusters_df.itertuples():
        try:
            # Initialize flag to track if all conditions passed
            passed_all_conditions = True
            edited = False
            
            # parse the EditingSites_to_PhredScore_Map and  MM_sites_map (presented as string)
            editing_sites_map = ast.literal_eval(read.EditingSites_to_PhredScore_Map)
            MM_sites_map = ast.literal_eval(read.MM_to_PhredScore_Map)
            total_num_of_ES = read.Number_of_Editing_Sites

            # Check conditions

            # filter ES and MM - count number of sites which qualities are bigger or equal to threshold. 
            # if the count is less than the minimum number of ES, the read did not passed this condition.
            # !!!!!!!!!!!!!!! ALL OTHER CONDITIONS WILL BE CHECKED ONLY ON THE PASSED ES & MM !!!!!!!!!!!!!!!!!!
            num_of_passed_ES, passed_ES_map = filter_sites(editing_sites_map, args.min_phred_score)
            num_of_passed_MM, passed_MM_map = filter_sites(MM_sites_map, args.min_phred_score)
            num_of_passed_MM += num_of_passed_ES
            min_phred_score_passed = True if num_of_passed_ES >= args.min_editing_sites else False

            # if no passed editing sites detected - continue to the next read
            # (keep column 'Edited' as False)
            if num_of_passed_ES == 0:
                continue
            # else - read was edited
            edited = True
            
            # cluster_len = last position of editing site list minus first position of editing sites list
            cluster_len = max(passed_ES_map.keys()) - min(passed_ES_map.keys())

            # editing fraction
            editing_fraction = 0 if num_of_passed_MM == 0 else num_of_passed_ES / num_of_passed_MM

            min_editing_sites_passed = check_Condition(num_of_passed_ES, args.min_editing_sites)
            min_editing_fraction_passed = check_Condition(editing_fraction, args.min_editing_fraction)
            min_es_length_ratio_passed = check_Condition(num_of_passed_ES / read.Alignment_length, args.min_es_length_ratio)
            min_cluster_length_ratio_passed = check_Condition(cluster_len / read.Alignment_length, args.min_cluster_length_ratio)

            #list of conditions -
            conditions_list = [min_editing_sites_passed, min_editing_fraction_passed, min_phred_score_passed, min_es_length_ratio_passed, min_cluster_length_ratio_passed]

            # Update passed_all_conditions flag
            passed_all_conditions = all(conditions_list)

            # if read passed all parameteres add a row to the passed_statisctics
            if passed_all_conditions and out_passed:
                avg_phred_score, avg_es_distance= es_statistics(passed_ES_map)
                es_statistics_data = {'Number_of_Total_ES': total_num_of_ES,'Number_of_Passed_ES': num_of_passed_ES, 'Number_of_Passed_MM': num_of_passed_MM,'Editing_Fraction_Passed':editing_fraction,'Passed_ES_Pos_to_Phred_Score_Map': passed_ES_map, 'Average_ES_Phred_Score': avg_phred_score, 'Average_Adjacent_ES_Distance': avg_es_distance}
                all_data = {**read._asdict(), **es_statistics_data} 
                passed_reads_data.append(all_data)
            # append details to the motifs file
            if passed_all_conditions and out_motifs:
                motifs_count_data = {"Read_ID": read.Read_ID, "Mate": read.Mate}
                motifs_count_data.update(motifs_count(read.Read_Sequence, passed_ES_map))
                motifs_data.append(motifs_count_data)

            if passed_all_conditions and out_bed:
                bed_rows = parse_bed_row(read, passed_ES_map)
                bed_file_data.extend(bed_rows)

            # append details to the condition_analysis file
            if out_analysis:
                condition_analysis_rows.append([read.Read_ID] + [passed_all_conditions] + [edited] + conditions_list)

        except Exception as e:
            print(f"Error processing read {read.Read_ID}: {e}")
            continue



def main():
    # read data from the detected clusters file
    clusters_df = pd.read_csv(sample_clusters_csv)
    multimappers_before_reads_count = clusters_df.shape[0]

    # deal with multimappers - 
    # 1- add editing fration based on Number_of_Editing_Sites & Number_of_total_MM
    # Step 1: Add editing fraction based on Number_of_Editing_Sites and Number_of_total_MM
    # Ensure the fraction is zero when Number_of_total_MM is zero
    clusters_df['Editing_Fraction'] = np.where(clusters_df['Number_of_total_MM'] == 0, 
                                                0, 
                                            clusters_df['Number_of_Editing_Sites'] / clusters_df['Number_of_total_MM'])
    # 2 - group using Read_ID and select the read with the highest fraction in the group 
    clusters_df = clusters_df.loc[clusters_df.groupby('Read_ID')['Editing_Fraction'].idxmax()]
    before_filter_reads_count = clusters_df.shape[0]

    if (before_filter_reads_count > 0):
        if (user_batch_size and (user_batch_size > 0)):
            batch_size = user_batch_size
        else:
            batch_size = int(before_filter_reads_count / max_threads)
        # Split the DataFrame into smaller batches
        batches = [clusters_df.iloc[i:i + batch_size] for i in range(0, len(clusters_df), batch_size)]

        # Use ThreadPoolExecutor to parallelize processing
        with ThreadPoolExecutor(max_workers=max_threads) as executor:  
            futures = [executor.submit(process_reads, batch) for batch in batches]
            for future in as_completed(futures):
                future.result()  # Ensure that exceptions are raised

    
    # WRITE OUTPUT
    #CREATE DATA FRAMES
    if(len(passed_reads_data) == 0):
        f = open(f'{sample_id}_no_passed_reads.txt', "w")
        sys.exit()
    

    # create passed DF and change order of columns
    passed_df = pd.DataFrame(passed_reads_data)
    column_order = [
    'Index', 'Read_ID', 'Mate', 'Chromosome', 'Strand', 'Position_0based', 'Alignment_length',
    'Read_Sequence', 'Visualize_Allignment', 'Reference_Sequence', 'cigar', 'flag',
    'Genomic_Position_Splicing_Blocks_0based', 'Read_Relative_Splicing_Blocks_0based',
    'Number_of_Total_ES', 'Number_of_Passed_ES', 'Number_of_total_MM', 'Number_of_Passed_MM',
    'Editing_Fraction_Passed', 'Passed_ES_Pos_to_Phred_Score_Map', 'Average_ES_Phred_Score',
    'Average_Adjacent_ES_Distance', 'A2C_MM', 'A2G_MM', 'A2T_MM', 'C2A_MM', 'C2G_MM',
    'C2T_MM', 'G2A_MM', 'G2C_MM', 'G2T_MM', 'T2A_MM', 'T2C_MM', 'T2G_MM', 'Ref2N_MM',
    'NtoAlt_MM'
    ]
    # Keep only the columns that exist in both the DataFrame and the specified order
    filtered_columns = [col for col in column_order if col in passed_df.columns]
    # Rearrange the columns in the DataFrame
    passed_df = passed_df[filtered_columns]

    # create condition df
    condition_header = [
        'Read_ID', 'Passed_All', 'Edited',
        'Min_Editing_Sites_' + str(args.min_editing_sites), 
        'Min_Editing_to_Total_MM_Fraction_' + str(args.min_editing_fraction), 
        'Min_Editing_Phred_Score_' + str(args.min_phred_score), 
        'Min_Editing_to_Read_Length_Ratio_' + str(args.min_es_length_ratio), 
        'Min_Cluster_Length_to_Read_Length_Ratio_' + str(args.min_cluster_length_ratio)
        ]
    condition_df = pd.DataFrame(condition_analysis_rows, columns=condition_header)
    
    motifs_df = pd.DataFrame(motifs_data)
    if out_passed:
        # Write all the rows to the passed CSV file
        passed_df.to_csv(passed_csv_path, index=False)
    if out_analysis:
        # Write all the rows to the condition_analysis CSV file
        condition_df.to_csv(condition_analysis_path, index=False)
    if out_motifs:
        # Write all the rows to the condition_analysis CSV file
        motifs_df.to_csv(motifs_csv_path, index=False)
    if out_bed:
        # Create a new DataFrame for the BED file
        bed_df = pd.DataFrame(bed_file_data)
        # Save to a BED file
        bed_df.to_csv(bed_path, sep='\t', header=False, index=False)
    
    if out_json:
        json_data = create_json(passed_df,motifs_df, condition_df,multimappers_before_reads_count, before_filter_reads_count)
        # write the data to a json file
        with open(json_summary_path, 'w') as f:
            json.dump(json_data, f, indent = 4)

if __name__ == "__main__":
    main()