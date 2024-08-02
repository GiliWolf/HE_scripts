#!/bin/bash

# Define the source directories and the target directory
target_base="/private10/Projects/Gili/HE_workdir/first_part/test_samples/GTEX_data/WholeBlood/new"

# List of filenames to copy
files=(
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-TSE9/GTEX-TSE9_WholeBlood_GTEX-TSE9-0005-SM-4DXUF_1.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-TSE9/GTEX-TSE9_WholeBlood_GTEX-TSE9-0005-SM-4DXUF_2.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-O5YT/GTEX-O5YT_WholeBlood_GTEX-O5YT-0007-SM-32PK7_1.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-O5YT/GTEX-O5YT_WholeBlood_GTEX-O5YT-0007-SM-32PK7_2.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-ZTX8/GTEX-ZTX8_WholeBlood_GTEX-ZTX8-0006-SM-4YCE4_1.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-ZTX8/GTEX-ZTX8_WholeBlood_GTEX-ZTX8-0006-SM-4YCE4_2.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-14ICK/GTEX-14ICK_WholeBlood_GTEX-14ICK-0006-SM-5NQB5_1.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-14ICK/GTEX-14ICK_WholeBlood_GTEX-14ICK-0006-SM-5NQB5_2.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-TMMY/GTEX-TMMY_WholeBlood_GTEX-TMMY-0005-SM-33HBN_1.fastq.bz2"
        "/private3/GTEx/RNA-Seq/WholeBlood/GTEX-TMMY/GTEX-TMMY_WholeBlood_GTEX-TMMY-0005-SM-33HBN_2.fastq.bz2"
)

# Copy each file from the source to the target directory
for file in "${files[@]}"; do
    cp "${file}" "${target_base}"
    echo "Finished copying $file"
    
    # Unzip the file in the target directory
    bunzip2 "${target_base}/$(basename ${file})"
    echo "Finished unzipping $(basename ${file})"
done

echo "Finished processing all files"

