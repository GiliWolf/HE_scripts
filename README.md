# RNA Hyper Editing Detection tool

<p> Hyper-edited RNA detection tool based on the algorithm presentes in:<br> A genome-wide map of hyper-edited RNA reveals numerous new sites, Hagit T. Porath, Shai Carmi & Erez Y. Levanon, Nature Communication 2014 (https://www.nature.com/articles/ncomms5726) </p>

___________________________________


## Dependecies:
  * Nextflow: https://github.com/nextflow-io/nextflow.git

___________________________________

## Main Parts:
  * [Genome Transformation](#Genome-Transformation)
  * [Main algorithm](#Main-algorithm)
  * [Detection](#Detection)

___________________________________

## Genome Transformation
<p> Transformation of the genome (12 transformations for each possible base combination) and 
indexing the transformed genomes (using STAR).</p>

### simple usage
>nextflow genome_setup.nf -c genome_setup.nf.config --genome_fasta <genome_fasta_file>

 - `-genome_fasta` : path to genome's fasta file to be transformed
- `-outdir` : path to the output directory

___________________________________
    
## Main algorithm
<p> Similar to the original pipeline, the process involves aligning the samples to the original genome, then transforming the unmapped reads and mapping them again to the transformed genome. Finally, the second mapped reads are retransformed to their original sequence.</p>

### simple usage
>nextflow -c pre_analysis.nf.config run pre_analysis.nf --reads_dir <reads_dir_path> --outdir <output_dir_path> [--pair_end]

### Input File Parameters

- `--pair_end`
  - Pair end flag (0: Single End, 1: Paired End)

- `--reads_dir`
  - Directory containing the samples

- `--reads_suffix`
  - Suffix of the sample files (default: ".fastq")

- `--file_seperator`
  - Separator of the file's name components (default: "_")

- `--mate_seperator`
  - Suffix of both mates (default: "_")

- `--mate1_suff`
  - Suffix of mate1 (default: "1")

- `--mate2_suff`
  - Suffix of mate2 (default: "2")

### Output Directories Parameters

- `--outdir`
  - Path for the output directory (relative path preferred)

- `--fastp_output_dir`
  - Path for the output of the fastp process (default: "${params.outdir}/fastp")

- `--first_map_output_dir`
  - Path for the output of the first map process (default: "${params.outdir}/first_map")

- `--transform_output_dir`
  - Path for the output of the transform reads process (default: "${params.outdir}/transformed_unmapped")

- `--second_map_output_dir`
  - Path for the output of the second map process (default: "${params.outdir}/second_map")

- `--retransform_output_dir`
  - Path for the output of the second map process (default: "${params.outdir}/re-transform")

### Fastp Parameters

- `--N_bases_num`
  - Maximum number of N bases in a read (default: 5)

- `--avg_quality`
  - Minimum average quality score of a read (default: 30)

- `--low_quality_per`
  - Minimum percentage of bases allowed to be unqualified (default: 20)

- `--low_quality_num`
  - Minimum quality value that a base is qualified, Phred score (default: 25)

### First Map Parameters

- `--fastp_output_suffix`
  - Suffix of the unmapped files (default: ".processed.fastq")

- `--read_files_command`
  - Command line to execute for each input file (default: "cat")

- `--genome_index_dir`
  - Path to the genome's index directory

- `--SAM_attr`
  - Desired SAM attributes (default: "All")

- `--outSAMtype`
  - Type of the output mapped files (default: "BAM Unsorted")

- `--min_SJ_overhang`
  - Minimum overhang (i.e., block size) for spliced alignments (default: 8)

- `--max_intron_size`
  - Maximum intron length (default value for STAR and HE: 1,000,000)

- `--max_mates_gap`
  - Maximum genomic distance between mates (default: 600,000)

- `--max_mismatches_ratio_to_ref`
  - Maximum ratio of mismatches to mapped length (default value for STAR and HE: 3)

- `--max_mismatche_ratio_to_read`
  - Maximum ratio of mismatches to read length (default value for STAR and HE: 1)

- `--norm_num_of_matches`
  - Minimum number of matched bases normalized to the read length (sum of mates’ lengths for PE reads) (default value for STAR and HE: 0.66)

- `--max_num_of_allignment`
  - Maximum number of multiple alignments allowed for a read; if exceeded, the read is unmapped (default: 1)

- `--genome_load_set`
  - Genome shared memory setting (default value for STAR and HE: "NoSharedMemory")

- `--num_of_threads`
  - Number of threads for each mapping process (default: 5)

- `--unmapped_out_files`
  - Unmapped reads will be output into separate file(s) `Unmapped.out.mate1[2]`, formatted the same way as input read files (default: "Fastx")

- `--output_files_permissions`
  - File permissions (default: "All_RWX")

### Second Map Parameters

- `--transformed_indexes`
  - Path to the directory of the transformed genome indexes from HE part 1 (e.g., "/generic_transform/hg38transform/genome_indexing/*")

- `--STAR_MAX_PARALLEL`
  - Number of files to be mapped in parallel (default: 6)

- `--second_map_genome_load_set`
  - Controls how the genome is loaded into memory (default: "LoadAndKeep")

### Retransform Parameters

- `--original_reads`
  - Path to the original fastq files (default: extracts from `first_map_output_dir`)

- `--filter_sam_files`
  - STAR output files to be output (default: '*Aligned.out*')

___________________________________

## Detection
1.  Analyze and detect the editing events from the overall general mismatches events.
	- CSV output file with:
		- *basic*: 'Read_ID', 'Chromosome', 'Strand', 'Position_0based','Alignment_length','Read_Sequence', 'Reference_Sequence','Number_of_MM', 'Number_of_Editing_Sites', 'Editing_to_Total_MM_Fraction' 
		- *all*: 'Read_ID', 'Chromosome', 'Strand', 'Position_0based','Alignment_length','Read_Sequence', 'Visualize_Allignment','Reference_Sequence', 'cigar', 'flag', 'Genomic_Position_Splicing_Blocks_0based','Read_Relative_Splicing_Blocks_0based', 'Number_of_total_MM', 'Number_of_Editing_Sites', 'Editing_to_Total_MM_Fraction', 'EditingSites_to_PhredScore_Map', 'MM_to_PhredScore_Map'
2. filter HE read based on pre-defined conditions
	- *Passed* CSV with the information of the reads that passed all of the consitions
	- *Condition analysis* CSV file with True/False for each of the conditions (passed/not passed) for each read
	- *Summary* JSON file of the sample

### simple usage
>./nextflow -c HE_detection.config.nf run HE_detection.n --reads_dir <path_to_transformed_sam_dir> --outdir <output_dir_path> [--pair_end]

>#for independent run:
>
>>./nextflow -c HE_detection.config.nf run HE_detection.nf --HE_reads <path_to_sam_file> --fasta_path <path_to_genome_fasta> --outdir <output_dir_path> [--pair_end] --ref_base <ref_base{A/G/C/T}> --alt_base <alt_base{A/G/C/T}> -entry independent
>
>>#for example: ./nextflow -c HE_detection.config.nf run HE_detection.nf --HE_reads /sam_dir/*.sam --fasta_path hg38.fa --outdir independent_HE --pair_end --ref_base A --alt_base G -entry independent
### Input paths

- ``reads_suffix``
  Suffix of the reads files.

- ``mate_seperator``
  Separator for mates in the file names.

- ``mate1_suff``
  Suffix for mate 1.

- ``mate2_suff``
  Suffix for mate 2.

- ``input_dir``
  Directory containing input files.

- ``SE_HE_reads``
  Path pattern for single-end reads.

- ``PE_HE_reads``
  Path pattern for paired-end reads.

- ``pair_end``
  Pair end flag (0: Single End, 1: Paired End).

- ``fasta_path``
  Path to the reference genome FASTA file.

- ``detect_python_script``
  Path to the Python script for detecting clusters.

- ``filter_python_script``
  Path to the Python script for filtering clusters.

- ``PE_filter_python_script``
  Path to the Python script for filtering PE clusters (if applicable).

- ``grid_serch_python_script``
  Path to the Python script for grid search.

- ``file_seperator``
  Separator in file names.

- ``python_command``
  Command to execute Python scripts.

### For Independent Run

- ``ref_base``
  Base reference for independent run.

- ``alt_base``
  Base alteration for independent run.

### Output Directories

- ``outdir``
  Base output directory path.

- ``detect_output_dir``
  Output directory path for detected clusters.

- ``filter_output_dir``
  Output directory path for filtered clusters.

- ``PE_filter_output_dir``
  Output directory path for PE filtered clusters.

- ``grid_search_output_dir``
  Output directory path for grid search results.

### Detect Script Arguments

- ``columns_select``
  Columns to select in detection script ('all' or specific columns).

### Filter PE

- ``unmapped_fastq``
  Path to the directory of unmapped fastq files.

- ``unmapped_star_file_format_mate1``
  Suffix format for unmapped mate 1 STAR files.

- ``unmapped_star_file_format_mate2``
  Suffix format for unmapped mate 2 STAR files.

### Ranges of Grid Search

#### Editing Sites

- ``es_start``
  Start value for editing sites.

- ``es_end``
  End value for editing sites.

- ``es_step``
  Step value for editing sites.

####  Editing Fraction

- ``ef_start``
  Start value for editing fraction.

- ``ef_end``
  End value for editing fraction.

- ``ef_step``
  Step value for editing fraction.

####  Phred Score

- ``ps_start``
  Start value for Phred score.

- ``ps_end``
  End value for Phred score.

- ``ps_step``
  Step value for Phred score.

####  Num of ES to Read Length Ratio

- ``es2l_start``
  Start value for the ratio of editing sites to read length.

- ``es2l_end``
  End value for the ratio of editing sites to read length.

- ``es2l_step``
  Step value for the ratio of editing sites to read length.

#### Cluster Length to Read Length Ratio

- ``cl2l_start``
  Start value for the ratio of cluster length to read length.

- ``cl2l_end``
  End value for the ratio of cluster length to read length.

- ``cl2l_step``
  Step value for the ratio of cluster length to read length.

### Merge JSON

- ``remove_jsons``
  Flag to remove JSONs.
