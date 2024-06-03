from collections import Counter
import pandas as pd

def sample_motif(nt_count_df, location):
    nt_precentage = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    num_reads, num_columns = nt_count_df.shape
    for nt in nt_precentage.keys():
        perc = (nt_count_df[f'{location}_{nt}'].sum() / num_reads)
        nt_precentage[nt] = perc
    return nt_precentage


motif_location = ['upstream', 'downstream']
def create_json(passed_df, condition_df):

    # init general data dictonary 
    json_data = {}

    # add total number of filtered reads
    # json_data["number of passed reads: "] = int(condition_df['Passed_All'].sum())

    # average es number
    average_num_of_ES = passed_df['number_of_ES'].mean()
    json_data["average number of ES: "] = average_num_of_ES

    # get sample ES motif
    for loc in motif_location:
        nt_count_df = passed_df[[f'{loc}_{nt}' for nt in ['A', 'C', 'G', 'T']]]
        nt_precentage =sample_motif(nt_count_df, loc)
        json_data[f'{loc}_motif: '] = nt_precentage

def es_statistics(editing_sites_map):
    id = "id1"
    phred_scores = editing_sites_map.values()
    phred_score_avg = sum(phred_scores) / len(phred_scores)
    es_positions = list(editing_sites_map.keys())
    es_positions.sort()
    avg_lambda_func =  lambda es_positions: sum(abs(es_positions[i+1] - es_positions[i]) for i in range(0, len(es_positions) - 1)) / len(es_positions)
    avg_es_distance = avg_lambda_func(es_positions)
    
    return phred_score_avg, avg_es_distance


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
reference_sequence="AACGGAAGGAA"
complementary = { 'A':'T', 'T':'A', 'G':'C','C':'G' }
rev_comp_seq = ''.join(([complementary[i] for i in reversed(reference_sequence)]))
print(rev_comp_seq)
map = {2: 36, 6: 36, 8: 36}
cluster_len = (15-6)
seq = "AACGGAAGGAA"

passed_read_data = []
all_data = {}
avg_phred_score, avg_es_distance, num_of_es = es_statistics(editing_sites_map)
motifs_count_data = motifs_count(read.Read_Sequence, editing_sites_map)
es_statistics_data = {'number_of_ES': num_of_es, 'average_es_phred_score': avg_phred_score, 'average_adjacent_es_distance': avg_es_distance}
all_data = {"Read_ID": "id"}
all_data.update(es_statistics_data)
all_data.update(motifs_count_data)
passed_read_data.append(all_data)


passed_df = pd.DataFrame(passed_read_data)
create_json(passed_df, 0)
