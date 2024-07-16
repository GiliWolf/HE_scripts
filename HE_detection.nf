// usage - ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf
// independent - 
// ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf --SE_HE_reads <path_to_sam> --fasta_path <path_to_genome_fasta> --outdir <outdir_path> --pair_end 0 --ref_base A --alt_base G -entry independent

process DETECT {
        maxForks 3
        tag "detection: ${sample_id}"
        publishDir "${params.detect_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(genome)
        path(python_script)

    output:
        tuple val(base_comb), val(sample_id), path('*detected.csv')
        // stdout

    script:
        outout_path = "${sample_id}${params.file_seperator}detected.csv"
        comb_bases_list = base_comb.tokenize("2")
        ref_base = comb_bases_list[0]
        alt_base = comb_bases_list[1]

        """
        ${params.python_command} ${python_script} -i ${file} -g ${genome} -o ${outout_path} -rb ${ref_base} -ab ${alt_base} -t ${params.max_detection_threads} -b ${params.detection_batch_size} -c ${params.detection_columns_select}

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
        analysis_output_path = "${sample_id}${params.file_seperator}condition_analysis.csv"
        summary_output_path = "${sample_id}${params.file_seperator}summary.json"
    """
    ${params.python_command} ${filter_python_script} -i ${file} -id ${sample_id} -o ${filtered_output_path} -O ${analysis_output_path} -j ${summary_output_path} -t ${params.max_filter_threads} -b ${params.filter_batch_size} -ot ${params.filter_output_types} -es ${params.min_editing_sites} -ef ${params.min_editing_fraction} -ps ${params.min_phred_score} -es2l ${params.min_es_length_ratio} -cl2l ${params.min_cluster_length_ratio}
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
 

process GRID_SEARCH_FILTER {
        maxForks 3
        tag "GS filter: ${sample_id}"
        publishDir "${params.grid_search_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)
        path(grid_search_python_script)

    output:
        tuple val(base_comb), path('*')
        
    script:
    """
    ${params.python_command} ${grid_search_python_script} -i ${file} -s ${filter_python_script} -bc ${base_comb} -id ${sample_id} -D "./" -es ${params.es_start} ${params.es_end} ${params.es_step} -ef ${params.ef_start} ${params.ef_end} ${params.ef_step} -ps ${params.ps_start} ${params.ps_end} ${params.ps_step} -es2l ${params.es2l_start} ${params.es2l_end} ${params.es2l_step} -cl2l ${params.cl2l_start} ${params.cl2l_end} ${params.cl2l_step}

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

workflow independent{

    def base_comb="${params.ref_base}2${params.alt_base}"

    if (params.pair_end == 0)                         //SE
        Channel
            .fromPath(params.SE_HE_reads, checkIfExists: true)
            .map {file -> tuple (base_comb,file.baseName, file)}
            .set {samples_ch}
    else if (params.pair_end == 1)                    //PE
        Channel
            .fromPath(params.PE_HE_reads, checkIfExists: true)
            .map {file -> tuple (base_comb, file.baseName, file)}
            .set {samples_ch} 
    else             // raise error if pair_end flag != 0/1
        error "----------------\n error: pair_end flag must be 0/1\n----------------"


    detecter_clusters_ch = DETECT(samples_ch, params.fasta_path, params.detect_python_script)
    
    if (!params.no_filter){
        after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script)
    }

    if (params.grid_search){
        grid_search_ch = GRID_SEARCH_FILTER(detecter_clusters_ch, params.GS_filter_script, params.grid_serch_python_script)
        merged_ch = MERGE_GRID_OUTPUT(grid_search_ch)
        files_to_remove = grid_search_ch.map {tuple -> tuple[1]}.filter (~/.*\.json$/)
        REMOVE_JSONS(files_to_remove)
    }
}




workflow {

    // GET SAMPLES:
    Channel
        .fromPath(params.HE_reads, checkIfExists: true)
        .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.baseName, file)}
        .set {samples_ch} 

    detecter_clusters_ch = DETECT(samples_ch, params.fasta_path, params.detect_python_script)
    
    if (!params.no_filter){
        after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script)
    }
    if (params.grid_search){
        grid_search_ch = GRID_SEARCH_FILTER(detecter_clusters_ch, params.GS_filter_script, params.grid_serch_python_script)
        merged_ch = MERGE_GRID_OUTPUT(grid_search_ch)
        files_to_remove = grid_search_ch.map {tuple -> tuple[1]}.filter (~/.*\.json$/)
        REMOVE_JSONS(files_to_remove)
    }
    // Channel
    //     .fromPath("/private10/Projects/Gili/HE_workdir/detection/BrainCerebellum_PE_multimappers/grid_search/**/*.json")
    //     .filter(~/.*(?<!merged\.json)$/)
    //     .set {files_to_remove_ch}

    // //files_to_remove_ch.view()
    // rm_ch = REMOVE_JSONS(files_to_remove_ch)
    // rm_ch.view()
}
