// usage - ./nextflow -c HE_scripts/HE_detection.config.nf run HE_scripts/HE_detection.nf
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

        //Usage: python detect_clusters.py <bam_path> <fasta_path> <output_path> <ref_base> <alt_base>
        """
        ${params.python_command} ${python_script} ${file} ${genome} ${outout_path} ${ref_base} ${alt_base}
        
        """ 
}

process FILTER {
        tag "filter: ${sample_id}"
        publishDir "${params.filter_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(base_comb), val(sample_id), path(file)
        path(filter_python_script)
        path(grid_serch_python_script)

    output:
        path('*')
        // stdout
    script:
        filtered_output_path = "${sample_id}${params.file_seperator}filtered.csv"
        analysis_output_path = "${sample_id}${params.file_seperator}condition_analysis.csv"
        // """
        // ${params.python_command} ${filter_python_script} -i ${file} -o ${filtered_output_path} -O ${analysis_output_path}

        
        // """  

        //      HE_grid_search.py [-h] -i INPUT_FILE -id SAMPLE_ID -O OUTPUT_DIR [-es X X X] [-ef X X X] [-ps X X X] [-es2l X X X] [-cl2l X X X]
    """
    ${params.python_command} ${grid_serch_python_script} -i ${file} -id ${sample_id} -D "./" \
        -es ${params.es_start} ${params.es_end} ${params.es_step} \
        -ef ${params.ef_start} ${params.ef_end} ${params.ef_step} \
        -ps ${params.ps_start} ${params.ps_end} ${params.ps_step} \
        -es2l ${params.es2l_start} ${params.es2l_end} ${params.es2l_step} \
        -cl2l ${params.cl2l_start} ${params.cl2l_end} ${params.cl2l_step}

    """
}

workflow {  
        // GET SAMPLES:
    // if (params.pair_end == 0)                         //SE
    //     Channel
    //         .fromPath(params.SE_HE_reads, checkIfExists: true)
    //         .map {file -> tuple ("A2C",file.baseName, file)}
    //         .set {samples_ch}
    // else if (params.pair_end == 1)                    //PE
    //     Channel
    //         .fromFilePairs(params.PE_HE_reads, checkIfExists: true)
    //         .set {samples_ch} 
    // else             // raise error if pair_end flag != 0/1
    //     error "----------------\n error: pair_end flag must be 0/1\n----------------"

    Channel
        .fromPath(params.SE_HE_reads, checkIfExists: true)
        .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.baseName, file)}
        // .collect(flat:false)
        .set {samples_ch}
    // samples_ch.view()
    detecter_clusters_ch = DETECT(samples_ch, params.fasta_path, params.detect_python_script)

    after_filter_ch = FILTER(detecter_clusters_ch, params.filter_python_script, params.grid_serch_python_script)

}