// usage - ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf
// independent - 
// ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf --SE_HE_reads <path_to_sam> --fasta_path <path_to_genome_fasta> --outdir <outdir_path> --pair_end 0 --ref_base A --alt_base G -entry independent

process DETECT {
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

        //detect_clusters.py [-h] -i BAM_PATH -g FASTA_PATH -o OUTPUT_PATH -rb REF_BASE -ab ALT_BASE [-c {all,basic}]
        """
        ${params.python_command} ${python_script} -i ${file} -g ${genome} -o ${outout_path} -rb ${ref_base} -ab ${alt_base} -c ${params.columns_select}

        
        """ 
}

process FILTER {
        tag "filter: ${sample_id}"
        publishDir "${params.filter_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)

    output:
        tuple val(base_comb), val(sample_id), path('*')
        // stdout
    script:
        filtered_output_path = "${sample_id}${params.file_seperator}filtered.csv"
        analysis_output_path = "${sample_id}${params.file_seperator}condition_analysis.csv"
        summary_output_path = "${sample_id}${params.file_seperator}summary.json"
    """
    ${params.python_command} ${filter_python_script} -i ${file} -id ${sample_id} -o ${filtered_output_path} -O ${analysis_output_path} -j ${summary_output_path} -t all
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
        tag "filter: ${sample_id}"
        publishDir "${params.grid_search_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)
        path(grid_search_python_script)

    output:
        path('*')
        // stdout
    script:
    """
    ${params.python_command} ${grid_search_python_script} -i ${file} -s ${filter_python_script} -bc ${base_comb} -id ${sample_id} -D "./" -es ${params.es_start} ${params.es_end} ${params.es_step} -ef ${params.ef_start} ${params.ef_end} ${params.ef_step} -ps ${params.ps_start} ${params.ps_end} ${params.ps_step} -es2l ${params.es2l_start} ${params.es2l_end} ${params.es2l_step} -cl2l ${params.cl2l_start} ${params.cl2l_end} ${params.cl2l_step}

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
    
    after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script)
    // // after_filter_ch = GRID_SEARCH_FILTER(detecter_clusters_ch, params.filter_python_script, params.grid_serch_python_script)

}


workflow {

    // GET SAMPLES:
    if (params.pair_end == 0)                         //SE
        Channel
            .fromPath(params.SE_HE_reads, checkIfExists: true)
            .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.baseName, file)}
            .set {samples_ch}
    else if (params.pair_end == 1)                    //PE
        Channel
            .fromPath(params.PE_HE_reads, checkIfExists: true)
            .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.baseName, file)}
            .set {samples_ch} 
    else             // raise error if pair_end flag != 0/1
        error "----------------\n error: pair_end flag must be 0/1\n----------------"


    detecter_clusters_ch = DETECT(samples_ch, params.fasta_path, params.detect_python_script)
    
    after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script)
    // // after_filter_ch = GRID_SEARCH_FILTER(detecter_clusters_ch, params.filter_python_script, params.grid_serch_python_script)

}