
log.info"""
        HYPER-EDITING PRE-ANALYSIS:
        ===========================
        1. preprocess reads using fastp
        2. 1st STAR map: map original reads using STAR and get unmapped reads
        3. 12 transformation of the *unmapped* reads 
        4. 2nd STAR map: maping tranformed reads to fitted transformed genome
        5. retransform the *mapped* reads
        """
bases=['A','G','C','T']

// FASTP preprocessing the reads, in order to remove lowquality reads
process FASTP {
    tag "preprocessing: ${sample_id}"
    publishDir "${params.fastp_output_dir}", pattern: '*.{fastq,json}', mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id),path('*.fastq')
        path('*.json')


    script:
        """
        # run fastp on each sample
            # SE (run on all files ending with the reads suffixs)
            if [ $params.pair_end = 0 ]; then
                    ${params.fastp_command} \
                    -n ${params.N_bases_num} \
                    -e ${params.avg_quality} \
                    -u ${params.low_quality_per} \
                    -q ${params.low_quality_num} \
                    -j "${sample_id}.fastp.json" \
                    --dont_eval_duplication \
                    --in1 "${reads}" \
                    -o "${sample_id}.processed.fastq"
            else
            #PE (run on all 1st mate files, and change the mate suffixs for 2nd mate)
                    ${params.fastp_command} \
                    -n ${params.N_bases_num} \
                    -e ${params.avg_quality} \
                    -u ${params.low_quality_per} \
                    -q ${params.low_quality_num} \
                    -j "${sample_id}.fastp.json" \
                    --dont_eval_duplication \
                    --in1 ${reads[0]} \
                    -o "${sample_id}_1.processed.fastq" \
                    --in2 ${reads[1]} \
                    -O "${sample_id}_2.processed.fastq"
            fi
        

        """

}

//  Mapping all the reads using STAR and the originial genome's index, 
// in order to get the *unmapped* reads to further processing
process FIRST_STAR_MAP {
 
    tag "1st MAP: ${sample_id}"
    publishDir "${params.first_map_output_dir}/${sample_id}", pattern: "*", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        path genome_index

    output:
        tuple val(sample_id), path('*Unmapped*')
        // path('*[*Unmapped*]')

    script:
        """
        ${params.STAR_command} --readFilesCommand ${params.read_files_command} --readFilesIn ${fastq[0]} ${fastq[1]} --genomeDir ${genome_index} --outSAMattributes ${params.SAM_attr} --outSAMtype ${params.outSAMtype} --alignSJoverhangMin ${params.min_SJ_overhang} --alignIntronMax ${params.max_intron_size} --alignMatesGapMax ${params.max_mates_gap} --outFilterMismatchNoverLmax ${params.max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${params.max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${params.norm_num_of_matches} --outFilterMultimapNmax ${params.max_num_of_allignment} --genomeLoad ${params.genome_load_set} --runThreadN ${params.num_of_threads} --outReadsUnmapped ${params.unmapped_out_files} --runDirPerm All_RWX
        """
}

// Transform the unmapped reads 12 times for each possible non-equale base-pairs (for exmp- A>C, G>T...)
process TRANSFORM_READS {
    tag "${ref_base} >> ${alt_base} transform on ${sample_id}"
    publishDir "${params.transform_output_dir}/${sample_id}/mate1/", pattern: '*_1.fastq', mode: 'copy'
    publishDir "${params.transform_output_dir}/${sample_id}/mate2/", pattern: '*_2.fastq', mode: 'copy'

    input:
        each ref_base
        each alt_base
        tuple val(sample_id), path(unmapped_fastq)

    output:
        tuple val("${ref_base}2${alt_base}"), path('*.fastq')
    
    when:
        ref_base != alt_base

    script:
        """
        # create lowercase ref and alt bases
            ref_base_lower=\$(echo "${ref_base}" | tr '[:upper:]' '[:lower:]')
            alt_base_lower=\$(echo "${alt_base}" | tr '[:upper:]' '[:lower:]')

        # create uppercase ref and alt bases
            ref_base_upper=\$(echo "${ref_base}" | tr '[:lower:]' '[:upper:]')
            alt_base_upper=\$(echo "${alt_base}" | tr '[:lower:]' '[:upper:]')
        
        # create output path for each mate
            transformed_output_path_mate1="${ref_base}2${alt_base}_${sample_id}_1.fastq"
            transformed_output_path_mate2="${ref_base}2${alt_base}_${sample_id}_2.fastq"
       
        # generic awk script for each mate
            awk_script='
                NR%4==2 {
                    gsub(ref_base_lower, alt_base_lower);
                    gsub(ref_base_upper, alt_base_upper)
                } 1
            '

        # transform 1st mate 
            awk -v ref_base_lower="\${ref_base_lower}" -v alt_base_lower="\${alt_base_lower}" \
            -v ref_base_upper="\${ref_base_upper}" -v alt_base_upper="\${alt_base_upper}" \
            "\${awk_script}" "${unmapped_fastq[0]}" > "\${transformed_output_path_mate1}"
            
        # transform 2nd mate
            awk -v ref_base_lower="\${ref_base_lower}" -v alt_base_lower="\${alt_base_lower}" \
            -v ref_base_upper="\${ref_base_upper}" -v alt_base_upper="\${alt_base_upper}" \
            "\${awk_script}" "${unmapped_fastq[1]}" > "\${transformed_output_path_mate2}"

        """
}

process SECOND_STAR_MAP {
    tag "2nd MAP: ${transformed_reads}"
    publishDir "${params.second_map_output_dir}", pattern: '*.{SAM, BAM}', mode: 'copy'

    input:
    tuple path(transformed_reads), path(transformed_genome)

    output:
    // path mapped_reads
    stdout

    script:
    """
    echo "transformed_reads: ${transformed_reads}"
    echo "transformed genome ${transformed_genome}"


    """
}

// process RETRANSFORM {
//     tag "re-transform: ${read}"
//     publishDir "${params.retransform_output_dir}", pattern: '*.{fq,fastq}', mode: 'copy'

//     input:
//     path read

//     output:
//     path re_transform_reads

//     script:
//     '''
//     '''
    
// }
bases_combination=[['A','C'],['A','G'],['A','T'],['C','A'],['C','G'],['C','T'],['G','A'],['G','C'],['G','T'],['T','A'],['T','C'],['T','G']]
// PE workflow
workflow {
//     Channel
//         .fromFilePairs(params.PE_reads, checkIfExists: true)
//         .set {read_pairs_ch} 
//     FASTP(read_pairs_ch)

//     Channel
//         .fromPath(params.genome_index_dir)
//         .set {basic_index_ch}
//     FIRST_STAR_MAP(FASTP.out[0], basic_index_ch)

//     trans_reads_ch = TRANSFORM_READS(bases, bases,FIRST_STAR_MAP.out[0])


//     // get all tramsformed genome and extract prefix - ref_base2alt_base, and map it to the file
//     Channel
//         .fromPath(params.transformed_genomes)
//         .map {file -> tuple(file.name.toString().tokenize('_').get(0), file)}
//         .set {transformed_genomes_ch}
    
//     // concat genome chanel and reads chanel and grupy by the base comination
//     // (for exp: [A2C, A2C_transformed_genome]+[A2C, A2C_read1, A2C,read2]) --> 
//     //           [A2C,[A2C_transformed_genome,[A2C, A2C_read1, A2C,read2]]])
//     temp_reads_and_genomes_ch = trans_reads_ch.concat(transformed_genomes_ch).groupTuple()

//     // get only transformed files
//     reads_and_genomes_ch = temp_reads_and_genomes_ch.map {it -> it[1]}
//     SECOND_STAR_MAP(reads_and_genomes_ch)
//     // SECOND_STAR_MAP.out.view()

//     Channel
//         .fromPath(params.transformed_indexes)
//         .set {index_ch}
//     index_ch.view()
    Channel
        .from(bases_combination)
        .set {bases_combination_ch}
    bases_combination_ch.view()

}

   // Channel
    //     // .fromFilePairs(params.fastq_reads, checkIfExists: true)
    //     // .set {fastq_pairs_ch}
    //     .fromFilePairs(params.fastq_reads, checkIfExists: true)
    //     .map { sample_id, files -> tuple(sample_id, files.toList()) }
    //     .collect()
    //     .set {fastq_pairs_ch}
    // fastq_pairs_ch.view()