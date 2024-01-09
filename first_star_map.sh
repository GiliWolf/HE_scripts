
STAR_command="STAR-2.7.3a"
read_files_command="cat" #command line to execute for each of the input file
reads="/home/alu/fulther/Validations/GiliTest/first_part/pre_analysis/fastp/SRR11548778_1.processed.fastq /home/alu/fulther/Validations/GiliTest/first_part/pre_analysis/fastp/SRR11548778_2.processed.fastq" #paths to files that contain input read1 [read2]
genome_index_dir="/private/dropbox/Genomes/Human/hg38.STAR.7.ReadsLn100.gencode28"
SAM_attr="All" #desired SAM attributes
outSAMtype="BAM Unsorted"
file_prefix="./SRR11548778/" #"${first_map_output_dir}${smaple_id}"
min_SJ_overhang=8 #minimum overhang (i.e. block size) for spliced alignments
max_intron_size=1000000 #maximum intron length (dafault)
max_mates_gap=600000 #maximum genomic distance between mates
max_mismatches_ratio_to_ref=0.3 #max ratio of mismatches to *mapped* length (dafualt)
max_mismatche_ratio_to_read=1 #max ratio of mismatches to *read* length (dafualt)
norm_num_of_matches=0.66 #min number of matched bases normalized to the read length (sum of mates’ lengths for PE reads) (default)
max_num_of_allignment=1 #max number of multiple alignments allowed for a read (if exceeded->unmapped)
genome_load_set="NoSharedMemory" #genome shared memory is not used (dafault)
num_of_threads=1 #?
unmapped_out_files="Fastx" # unmapped reads will be output into separate file(s)Unmapped.out.mate1[2], formatted the same way as input read files

${STAR_command} --readFilesCommand ${read_files_command} --readFilesIn ${reads} --genomeDir ${genome_index_dir} --outSAMattributes ${SAM_attr} --outSAMtype ${outSAMtype} --outFileNamePrefix ${file_prefix} --alignSJoverhangMin ${min_SJ_overhang} --alignIntronMax ${max_intron_size} --alignMatesGapMax ${max_mates_gap} --outFilterMismatchNoverLmax ${max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${norm_num_of_matches} --outFilterMultimapNmax ${max_num_of_allignment} --genomeLoad ${genome_load_set} --runThreadN ${num_of_threads} --outReadsUnmapped ${unmapped_out_files}

