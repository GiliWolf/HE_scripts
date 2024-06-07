#!/bin/bash

# Define the source directories and the target directory
source_base="/private9/Projects/dsRNAProject/HE/ArteryAorta/unMap/"
target_base="/private10/Projects/Gili/HE_workdir/first_part/test_samples/GTEX_data/ArteryAorta"

# List of directories to process
dirs=("BrainCerebellum" "MuscleSkeletal" "WholeBlood")

# List of filenames to copy
files=(
    "GTEX-RUSQ_ArteryAorta_GTEX-RUSQ-0326-SM-47JWS-1.mem.um.fastq"
    "GTEX-RUSQ_ArteryAorta_GTEX-RUSQ-0326-SM-47JWS-2.mem.um.fastq"
    "GTEX-T6MN_ArteryAorta_GTEX-T6MN-1126-SM-4DM71-1.mem.um.fastq"
    "GTEX-T6MN_ArteryAorta_GTEX-T6MN-1126-SM-4DM71-2.mem.um.fastq"
    "GTEX-V955_ArteryAorta_GTEX-V955-0926-SM-4JBJ8-1.mem.um.fastq"
    "GTEX-V955_ArteryAorta_GTEX-V955-0926-SM-4JBJ8-2.mem.um.fastq"
    "GTEX-XUW1_ArteryAorta_GTEX-XUW1-1126-SM-4BONZ-1.mem.um.fastq"
    "GTEX-XUW1_ArteryAorta_GTEX-XUW1-1126-SM-4BONZ-2.mem.um.fastq"
    "GTEX-XV7Q_ArteryAorta_GTEX-XV7Q-0526-SM-4BRWR-1.mem.um.fastq"
    "GTEX-XV7Q_ArteryAorta_GTEX-XV7Q-0526-SM-4BRWR-2.mem.um.fastq"
)

# Copy each file from the source to the target directory
for file in "${files[@]}"; do
    cp "${source_base}${file}" "${target_base}"
    echo "finished coping $file"
done
echo "finished"


# # Loop through each directory
# for dir in "${dirs[@]}"; do
#   # Create target directory if it doesn't exist
#   target_dir="$target_base/$dir"
#   mkdir -p "$target_dir"
#   # # Loop through each file in the directory
#   #   for file in "$target_dir"/*; do
#   #   # Check if it is a file (not a directory)
#   #   if [[ -f "$file" ]]; then
#   #       decompressed_file="$target_dir/$(basename "${file%.bz2}")"
#   #       bzcat "$file" > "$decompressed_file"
#   #   fi
#   #   done
#   # Find the first 10 files and create symbolic links in the target directory
#   find "$source_base/$dir" -maxdepth 2 -type f | head -n 10 | while read -r file; do
#     #ln -s "$file" "$target_dir"
#     if [[ "$file" == *.bz2 ]]; then
#       decompressed_file="$target_dir/$(basename "${file%.bz2}")"
#       bzcat "$file" > "$decompressed_file"
#       echo "finished decompressing $decompressed_file"
#     fi
#   done
# done
