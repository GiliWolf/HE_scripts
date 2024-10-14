/*
----------------------------
Author: Gili Wolf
Date: 01-08-2024
Affiliation: Levanon Lab
----------------------------

Description: 
    This script is the third part of the Hyper Editing tool, inspired by Hagit T. Porath’s 2014 work ("A genome-wide map of hyper-edited RNA reveals numerous new sites", Nature Communications). The tool is designed to identify hyper-edited A-to-G clusters.

    First, the script processes BAM and BAI (indexed BAM) files to detect editing events by analyzing mismatches between read sequences and the reference genome, outputting all findings in the detect CSV file. Next, filtering is applied to identify hyper-edited reads based on Phred scores of editing events and mismatches, as well as criteria like editing fraction and cluster length. The filtering process generates several output files.

    **Important Notes:**
    - Multimappers: In cases of multiple loci, the locus with the highest editing scores (before Phred score filtering) is selected.
    - Conditions: All conditions are applied after Phred score filtering. Only editing sites (ES) and mismatches (MM) that pass the minimum Phred score are considered.

Usage: 
    ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf --detect_input_dir <DETECT_INPUT_PATH|ALIGN_OUTPUT_DIR_PATH> --detect_outdir <DETECT_OUTPUT_PATH> --genome_fasta <PATH_TO_GENOME> 
    --genome_index_dir <ORIGINAL_GENOME_INDEX_PATH>

    independent run - use --independent flag (and --index if BAI files are missing)

output files:
    - Condition Analysis CSV: This file provides a True/False assessment for each condition (passed/not
    passed) applied to each read.
    - Passed Reads CSV: Contains all the metadata provided in the initial detection CSV file, along
    with additional information post-filtering, such as the number of passed editing sites (ES) and
    mismatches (MM).
    - Motifs CSV: For each read, this file includes the count of each upstream base present in the
    total editing events (e.g., 3 of the events had T before them) and similarly for downstream
    bases. This information is used for later analysis of motifs related to the editing.
    - Clusters BED: Contains the genomic positions of each hyper-editing cluster.
    - Summary JSON: Provides statistics and summary information for the entire sample, such as
    the total number of passed reads and the average number of editing sites (ES)

Dependencies (using Docker containers):
   - samtools
   - python
   - PYSAM python module (version: 0.22.0)

----------------------------


*/

params.help=false

def helpMessage() {
    log.info '''
        HYPER-EDITING DETECTION:
        ===========================
        1. Detect all editing sites and mismatches using reference genome.
        2. Filter reads according to several criteria regarding ES and MM.
        3. Output statistics and summary files.
        ===========================
        usage: nextflow -c HE_detection.nf.config [options...] HE_detection.nf
        * for independent run (outsourced BAM files) :
          nextflow -c HE_detection.nf.config [options...] HE_detection.nf --independent
        * for help message:
          nextflow -c run HE_detection.nf --help

        options:
        ---------------
                --help              displays help message             
        ---------------
        Input File Parameters
                --bam_suffix        Suffix for input BAM files (default: ".bam")
                --bai_suffix        Suffix for input BAI (BAM index) files (default: ".bai")
                --detect_input_dir  Directory for detection input. Can be either the output directory from Pre-Analysis or a directory containing input BAM and/or BAI files (required)
                --retransform_dir_name  Directory where re-transformed BAM files from Pre-Analysis are stored (default: "re-transform"). Not used in an independent run
                --HE_reads          Path pattern for hyper-edited (HE) reads for a regular run (default: "$params.detect_input_dir/$params.retransform_dir_name/**/*$params.bam_suffix,$params.bai_suffix")
                --HE_reads_independent  Path pattern for HE reads in an independent run (default: "$params.detect_input_dir/**{$params.bam_suffix,$params.bai_suffix}")
                --detect_python_script  Path to the Python script for parallel detection (default: "/private10/Projects/Gili/HEworkdir/HEscripts/paralleldetection.py")
                --filter_python_script  Path to the Python script for parallel filtering (default: "/private10/Projects/Gili/HEworkdir/HEscripts/parallelfilter.py")
                --PE_filter_python_script  Path to the Python script for paired-end filtering
        ---------------
        Independent Run Options
                --independent       Boolean flag indicating whether the run is independent (default: false)
                --index             Boolean flag indicating whether indexing is enabled (default: false)
                --index_threads     Number of threads for indexing (default: 16)
        ---------------
        Output Parameters
                --detect_out_dir    Directory for storing detection outputs (required)
                --index_output_dir  Directory for storing index files (default: "$params.detect_out_dir/index")
                --detect_output_dir Directory for storing detected clusters files (default: "$params.detect_out_dir/detectedclusters")
                --filter_output_dir Directory for storing filtered clusters files (default: "$params.detect_out_dir/filteredclusters")
        ---------------
        Detection Script Parameters
                --detection_columns_select  Columns to include in detection output, either 'all' or 'basic' (default: 'all')
                --max_detection_threads  Maximum threads for parallel detection (default: 3)
                --detection_batch_size    Batch size for detection; set to 0 for automatic determination (default: 0)
        ---------------
        Filter Script Parameters
                --filter_output_types  Types of output for filtering: "all", "passed", "analysis", "motifs", "bed", or "summary" (default: "all")
                --max_filter_threads  Maximum threads for parallel filtering (default: 3)
                --filter_batch_size    Batch size for filtering; set to 0 for automatic determination (default: 0)
                --min_editing_sites    Minimum editing sites required for a read to be considered edited (default: 1)
                --min_editing_fraction Minimum fraction of editing sites relative to total sites for hyper-editing (default: 0.6)
                --min_phred_es_score   Minimum Phred score required for a base to be considered (default: 30)
                --length_ratio         Minimum ratio of editing site length to read length (default: 0.05)
                --min_cluster_length_ratio  Minimum ratio of cluster length to total read length (default: 0.1)

    '''
    .stripIndent()
}

process INDEX_BAM {
        maxForks 1
        tag "index: ${sample_id}"
        publishDir "${params.index_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(bam_file)

    output:
        tuple val(base_comb), val(sample_id), path(sorted_bam_path), path("*.bai"), env(COUNT)

    script:
        sorted_bam_path = "${sample_id}_sorted.bam"

        """
        samtools sort -o ${sorted_bam_path} ${bam_file}
        samtools index -@ ${params.index_threads} ${sorted_bam_path}
        COUNT=\$(samtools view -c ${sorted_bam_path})
        """ 
}

process COUNT_RECORDS {
        tag "count_records: ${sample_id}"
    input:
        tuple val(base_comb), val(sample_id), path(bam_file), path(index_bam_file)
    output:
        tuple val(base_comb), val(sample_id), path(bam_file), path(index_bam_file), env(COUNT)
    
    script:

        """
        COUNT=\$(samtools view -c ${bam_file})
        """ 
}
process DETECT {
        maxForks 3
        tag "detection: ${base_comb} -  ${sample_id}"
        publishDir "${params.detect_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(bam_file), path(bam_index_file), val(reads_count)
        path(genome)
        path(genome_index)
        path(python_script)

    output:
        tuple val(base_comb), val(sample_id), path('*detected.csv')

    script:
        outout_path = "${sample_id}${params.file_seperator}detected.csv"
        comb_bases_list = base_comb.tokenize("2")
        ref_base = comb_bases_list[0]
        alt_base = comb_bases_list[1]

        """
        ${params.python_command} ${python_script} -i ${bam_file} -f ${genome} -I ${bam_index_file} -F ${genome_index} -o ${outout_path} -rb ${ref_base} -ab ${alt_base} -t ${params.max_detection_threads} -n ${reads_count} -b ${params.detection_batch_size} -c ${params.detection_columns_select}

        """ 
}

process FILTER {
        maxForks 3
        tag "filter: ${sample_id}"
        publishDir "${params.filter_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)

    output:
        tuple val(base_comb), val(sample_id), path('*')
        // stdout
    script:
        filtered_output_path = "${sample_id}${params.file_seperator}passed.csv"
        motifs_output_path = "${sample_id}${params.file_seperator}motifs.csv"
        bed_output_path = "${sample_id}${params.file_seperator}clusters.bed"
        analysis_output_path = "${sample_id}${params.file_seperator}condition_analysis.csv"
        summary_output_path = "${sample_id}${params.file_seperator}summary.json"
    """
    ${params.python_command} ${filter_python_script} -i ${file} -id ${sample_id} -o ${filtered_output_path} -O ${analysis_output_path} -M ${motifs_output_path} -B ${bed_output_path} -j ${summary_output_path} -t ${params.max_filter_threads} -b ${params.filter_batch_size} -ot ${params.filter_output_types} -es ${params.min_editing_sites} -ef ${params.min_editing_fraction} -ps ${params.min_phred_score} -es2l ${params.min_es_length_ratio} -cl2l ${params.min_cluster_length_ratio}
    """
}

process PE_FILTER {
    tag "PE_filter: ${sample_id}"
    publishDir "${params.PE_filter_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(PE_python_script)

    output:
        path('*')
    
    script:
    """
    """
}
 

process GRID_SEARCH{
        maxForks 3
        tag "GS filter: ${sample_id}"
        publishDir "${params.grid_search_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)
        path(grid_search_python_script)

    output:
        tuple val(base_comb), path('*.json'), optional: true

    script:
    """
    ${params.python_command} ${grid_search_python_script} -i ${file} -s ${filter_python_script} -bc ${base_comb} -id ${sample_id} -th ${params.GS_max_threads} -D "" -es ${params.es_start} ${params.es_end} ${params.es_step} -ef ${params.ef_start} ${params.ef_end} ${params.ef_step} -ps ${params.ps_start} ${params.ps_end} ${params.ps_step} -es2l ${params.es2l_start} ${params.es2l_end} ${params.es2l_step} -cl2l ${params.cl2l_start} ${params.cl2l_end} ${params.cl2l_step}

    """
}

process MERGE_GRID_OUTPUT{
        tag "Merge jsons: ${base_comb}"
        publishDir "${params.grid_search_output_dir}/${base_comb}", pattern: '*_merged.json', mode: 'copy'

    input:
        tuple val(base_comb), path(json_files)

    output:
        path('*')
        stdout
    script:
    merged_output_file = "${base_comb}_merged.json"

    """
    # Check for the -remove flag
    REMOVE=false
    if [ ${params.remove_jsons} -eq 1 ]; then
        REMOVE=true
    fi

    # MODIFY FASTQS FILES INTO A BASH LIST:  
        # transform fastq file's names into string
        jsons_string="${json_files}"
        # Split the string into a list
        IFS=' ' read -ra json_list <<< "\$jsons_string"

    # Define the output file
    OUTPUT_FILE="${merged_output_file}"

    # Create or empty the output file
    > "\$OUTPUT_FILE"

    # Start the JSON array
    echo "[" > "\$OUTPUT_FILE"

    # Merge the JSON files
    for FILE in "\${json_list[@]}"; do
        echo "merging file: \$FILE"
        cat "\$FILE" >> "\$OUTPUT_FILE"
        echo "," >> "\$OUTPUT_FILE"
    done

    # Remove the trailing comma and add the closing bracket
    truncate -s-2 "\$OUTPUT_FILE"
    echo "]" >> "\$OUTPUT_FILE"

    # Remove individual JSON files if the -remove flag is set
    if [ "\$REMOVE" = true ]; then
        for FILE in "\${json_list[@]}"; do
            echo "removing file \$FILE:"
            rm \$FILE
        done
    fi

    """
}

process REMOVE_JSONS {
        tag "remove jsons: "

    input:
        path(json_files)

    output:
        stdout
    
    script:
    """
    # MODIFY FASTQS FILES INTO A BASH LIST:  
    # transform fastq file's names into string
    jsons_string="${json_files}"
    # Split the string into a list
    IFS=' ' read -ra json_list <<< "\$jsons_string"
    
    for FILE in "\${json_list[@]}"; do
        echo "removing file \$FILE:"
        rm \$(readlink -f \$FILE)
    done
    """

}


workflow HE_DETECTION {
    take:
        samples_ch
    
    main:
        if (params.index){
            index_ch = INDEX_BAM(samples_ch) 
            detected_clusters_ch = DETECT(index_ch, params.genome_fasta, params.genome_index_dir, params.detect_python_script)
        }
        else {
            samples_ch = COUNT_RECORDS(samples_ch)
            detected_clusters_ch = DETECT(samples_ch, params.genome_fasta, params.genome_index_dir, params.detect_python_script)
        }

        if (!params.no_filter){
            after_filter_ch = FILTER(detected_clusters_ch, params.filter_python_script)
        }
        if (params.grid_search){
            grid_search_ch = GRID_SEARCH(detected_clusters_ch, params.filter_python_script, params.grid_search_python_script)
            grid_search_ch.view()
            merged_ch = MERGE_GRID_OUTPUT(grid_search_ch)
            files_to_remove = grid_search_ch.map {tuple -> tuple[1]}
            REMOVE_JSONS(files_to_remove)
        }
}



workflow {
    if(params.help){
        helpMessage()
        System.exit()
    }

    if (params.independent){
        def base_comb="${params.ref_base}2${params.alt_base}"
        Channel
            .fromList(['A2C', 'A2G', 'A2T', 'C2A', 'C2G', 'C2T', 'G2A', 'G2C', 'G2T', 'T2A', 'T2C', 'T2G'])
            .set {bases_combination_ch}

        Channel
            .fromPath(params.HE_reads_independent, checkIfExists: true)
            .map {file -> tuple (file.name.toString().tokenize(params.suffix_seperator).get(0), file)}
            .groupTuple(by: 0, sort:true)
            .map {tuple -> tuple.flatten()}
            .set {files_ch}
        
        samples_ch = bases_combination_ch.combine(files_ch)
        samples_ch.view()
    }
    else {
        Channel
            .fromPath(params.HE_reads, checkIfExists: true)
            .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.name.toString().tokenize(params.suffix_seperator).get(0), file)}
            .groupTuple(by: [0,1], sort:true)
            .map {tuple -> tuple.flatten()}
            .set {samples_ch}
    }
    HE_DETECTION(samples_ch)
}