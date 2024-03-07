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
        // def ref_base = comb_bases_list[0]
        // def alt_base = comb_bases_list[1]
        // def comb_bases_list = ['A', 'C']
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
        path(python_script)

    output:
        path('*')
        // stdout
    script:
        filtered_output_path = "${sample_id}${params.file_seperator}filtered.csv"
        rejected_output_path = "${sample_id}${params.file_seperator}rejected.csv"
        """
        ${params.python_command} ${python_script} ${file} ${filtered_output_path} ${rejected_output_path}

        
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
    after_filtiration_ch = FILTER(detecter_clusters_ch, params.filter_python_script)

}