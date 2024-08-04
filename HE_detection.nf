// usage - ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf
// independent - 
// ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf --SE_HE_reads <path_to_sam> --fasta_path <path_to_genome_fasta> --outdir <outdir_path> --pair_end 0 --ref_base A --alt_base G -entry independent
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
        bed_output_path = "${sample_id}${params.file_seperator}es.bed"
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

// workflow independent{

//     def base_comb="${params.ref_base}2${params.alt_base}"

//     Channel
//             .fromPath(params.HE_reads, checkIfExists: true)
//             .map {file -> tuple (base_comb,file.baseName, file)}
//             .set {samples_ch}

//     index_ch = INDEX_BAM(samples_ch)

//     detecter_clusters_ch = DETECT(index_ch, params.fasta_path, params.detect_python_script)
    
//     if (!params.no_filter){
//         after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script)
//     }

//     if (params.grid_search){
//         grid_search_ch = GRID_SEARCH(detecter_clusters_ch, params.GS_filter_script, params.grid_search_python_script)
//         merged_ch = MERGE_GRID_OUTPUT(grid_search_ch)
//         files_to_remove = grid_search_ch.map {tuple -> tuple[1]}.filter (~/.*\.json$/)
//         REMOVE_JSONS(files_to_remove)
//     }
// }

workflow HE_DETECTION {
    take:
        samples_ch
    
    main: 
        index_ch = INDEX_BAM(samples_ch) 

        detected_clusters_ch = DETECT(index_ch, params.fasta_path, params.fasta_index_path, params.detect_python_script)
        
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

    if (params.independent){
        def base_comb="${params.ref_base}2${params.alt_base}"
        Channel
            .fromList(['A2C', 'A2G', 'A2T', 'C2A', 'C2G', 'C2T', 'G2A', 'G2C', 'G2T', 'T2A', 'T2C', 'T2G'])
            .set {bases_combination_ch}

        Channel
            .fromPath(params.HE_reads, checkIfExists: true)
            .map {file -> tuple (file.baseName, file)}
            .set {files_ch}
        
        samples_ch = bases_combination_ch.combine(files_ch)
    }
    else {
        Channel
            .fromPath(params.HE_reads, checkIfExists: true)
            .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.baseName, file)}
            .set {samples_ch}
    }

    HE_DETECTION(samples_ch)
}