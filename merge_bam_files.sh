#!/bin/bash
# nohup ./HE_scripts/merge_bam_files.sh &>> merge_ArteryAorta.out &
# Define variables
tissue="ArteryAorta"
output_dir="/private10/Projects/Gili/HE_workdir/HE_compare/merge_bam/${tissue}"

#mkdir $output_dir

# Define the list of sample IDs
samples=(
  "GTEX-111YS_ArteryAorta_GTEX-111YS-0526-SM-5GZXJ"
  "GTEX-139UW_ArteryAorta_GTEX-139UW-0426-SM-5K7V4"
  "GTEX-1212Z_ArteryAorta_GTEX-1212Z-0826-SM-5EQ51"
  "GTEX-O5YT_ArteryAorta_GTEX-O5YT-0426-SM-3MJHD"
  "GTEX-Q2AH_ArteryAorta_GTEX-Q2AH-0326-SM-48U1K"
  "GTEX-RUSQ_ArteryAorta_GTEX-RUSQ-0326-SM-47JWS"
  "GTEX-T6MN_ArteryAorta_GTEX-T6MN-1126-SM-4DM71"
  "GTEX-V955_ArteryAorta_GTEX-V955-0926-SM-4JBJ8"
  "GTEX-XUW1_ArteryAorta_GTEX-XUW1-1126-SM-4BONZ"
  "GTEX-XV7Q_ArteryAorta_GTEX-XV7Q-0526-SM-4BRWR"
)

echo $'////////////////////\n merging: ${tissue} \n ///////////////////////////'
# Loop over each sample ID
for sample_id in "${samples[@]}"; do
    first_map_bam="/private10/Projects/Gili/HE_workdir/first_part/GTEX_Multimappers_Reports/${tissue}/first_map/${sample_id}.Aligned.out.bam"
    re_transformed_bam="/private10/Projects/Gili/HE_workdir/first_part/GTEX_Multimappers_Reports/${tissue}/re-transform/A2G/A2G_${sample_id}_re-transformed.bam"
    
    # Define the output merged BAM file
    merged_bam="${output_dir}/${sample_id}_merged.bam"
  
    # Merge the BAM files using samtools
    echo "Merging BAM files for sample: $sample_id"
    samtools merge -f "$merged_bam" "$first_map_bam" "$re_transformed_bam"
  
    # Index the merged BAM file
    echo "Indexing merged BAM: $merged_bam"
    samtools index -@ 16 "$merged_bam"
done

echo $'////////////////////\n finished merging: ${tissue} \n ///////////////////////////'