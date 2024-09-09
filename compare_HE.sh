#!/bin/bash

# Define variables
tissue="WholeBlood"
original_tool_dir="/private9/Projects/dsRNAProject/HE/${tissue}/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2"
new_tool_dir="/private10/Projects/Gili/HE_workdir/detection/GTEX_Multimappers_All/${tissue}/filtered_clusters"
output_dir="/private10/Projects/Gili/HE_workdir/HE_compare/both/${tissue}"
#/private9/Projects/dsRNAProject/HE/ArteryAorta/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/GTEX-111YS_ArteryAorta_GTEX-111YS-0526-SM-5GZXJ.ES.bed_files/A2G.bed

# """
# for each sub directory:
# 1 - extract name and save it as "base_comb"
# 2 - get all files in the sub dir that ends with '_es.bed'
# for each file:
#     a - extract sample id by removing "[A-Z]2[A-Z]_" prefix and "_re-transformed_es.bed" suffix
#     b - find the original file in the original_tool_dir by the pattern: $original_tool_dir/$sample_ID.UE.bed_files/$base_comb.bed
#     c - run : bedtools intersect -a $file -b $original_file -v > ${base_comb}_${sample_id_compare}.bed
# """

mkdir -p "${output_dir}"

# Loop through each subdirectory in the new_tool_dir
for sub_dir in ${new_tool_dir}/*/; do
    # Extract the base directory name
    base_comb=$(basename "${sub_dir}")
    
    echo "base_comb: $base_comb"
    # Find all files ending with 'passed.csv' in the subdirectory
    for file in ${sub_dir}/*.bed; do
        # Extract sample ID by removing the specified prefixes and suffixes
        sample_id=$(basename "${file}" | sed -E 's/^[A-Z]2[A-Z]_//; s/_re-transformed_es\.bed$//')
        echo "sample_id: $sample_id"
        # Find the original file in the original_tool_dir
        original_file="${original_tool_dir}/${sample_id}.UE.bed_files/${base_comb}.bed"
        
        # Ensure the original file exists
        if [ -f "${original_file}" ]; then
            # Define the output file name
            output_file="${output_dir}/${base_comb}_${sample_id}_compare.bed"

            touch "${output_file}"
            # Run bedtools intersect
            bedtools intersect -a "${file}" -b "${original_file}" -wa > "${output_file}"
            
            echo "Processed ${file} and ${original_file}, output saved to ${output_file}"
        else
            echo "Original file ${original_file} not found, skipping."
        fi
    done
done