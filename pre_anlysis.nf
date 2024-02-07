/*
----------------------------
Author: Gili Wolf
Date: 01-02-2024
Affiliation: Levanon Lab
----------------------------
Description:
    This script is the second part of the Hyper Editing tool, inspired by the work of Hagit T. Porath in 2014 ("A genome-wide map of hyper-edited RNA reveals numerous new sites", Nature Communications). The aim of this tool is to identify hyper-edited A2G clusters.

    Similar to the original pipeline, the process involves aligning the samples to the original genome, then transforming the unmapped reads and mapping them again to the transformed genome. Finally, the second mapped reads are retransformed to their original sequence.

Usage: 
    nextflow -c PWD/pre_anylis.nf.config run PWD/pre_anlysis.nf

Dependencies (using Docker containers):
   - STAR aligner (version 2.7.10b)
   - FASTP preprocessor (version 0.23.4)
   - PYSAM python module (version: 0.22.0)

 Hardware Requirements:
   - Minimum RAM: ? GB
   - Minimum CPU: ? cores

 Additional Notes:
   - 

Disclaimer: This script is provided as-is without any warranty. Use at your own risk.
----------------------------

*/

log.info'''
        HYPER-EDITING PRE-ANALYSIS:
        ===========================
        1. preprocess reads using fastp
        2. 1st STAR map: map original reads using STAR and get unmapped reads
        3. 12 transformation of the *unmapped* reads 
        4. 2nd STAR map: maping tranformed reads to fitted transformed genome
        5. retransform the *mapped* reads
        '''
bases=['A','G','C','T']


/*
    Preprocess the reads using FASTP to filter out low-quality reads.
*/
process FASTP{
    tag "preprocessing: ${sample_id}"
    publishDir "${params.fastp_output_dir}", pattern: '*.{fastq,json}', mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path('*.fastq')
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

/*
    Map all the reads using STAR with the original genome's index to identify the unmapped reads for further processing.
*/
process FIRST_STAR_MAP{
    maxForks 1
    tag "1st MAP: ${genome_index}"
    publishDir "${params.first_map_output_dir}/${sample_id}", pattern: "*", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        path(genome_index)

    output:
        tuple val(sample_id), path('*Unmapped*')

    script:
        """
        ${params.STAR_command} --readFilesCommand ${params.read_files_command} --readFilesIn ${fastq[0]} ${fastq[1]} --genomeDir ${genome_index} --outSAMattributes ${params.SAM_attr} --outSAMtype ${params.outSAMtype} --alignSJoverhangMin ${params.min_SJ_overhang} --alignIntronMax ${params.max_intron_size} --alignMatesGapMax ${params.max_mates_gap} --outFilterMismatchNoverLmax ${params.max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${params.max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${params.norm_num_of_matches} --outFilterMultimapNmax ${params.max_num_of_allignment} --genomeLoad ${params.genome_load_set} --runThreadN ${params.num_of_threads} --outReadsUnmapped ${params.unmapped_out_files} --runDirPerm All_RWX
        """
}

/*
    Transform the unmapped reads 12 times for each possible non-equal base-pair substitution (e.g., A>C, G>T, etc.).
*/
process TRANSFORM_READS {
    tag "${ref_base} >> ${alt_base} transform on ${sample_id}"
    publishDir "${params.transform_output_dir}/${ref_base}2${alt_base}", pattern: '*.fastq', mode: 'copy'

    input:
        tuple val(sample_id), path(unmapped_fastq)
        each bases_combination

    output:
        path('*.fastq')
    
    script:
        ref_base = bases_combination[0]
        alt_base = bases_combination[1]

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

/*
    This process operates on one base combination (MM) at a time (A2C, C2A, T2G, ...).
    It takes in the transformed index of a certain base combination and all the transformed fastq files of the same combination.
    The process then maps each sample to the index using internal bash parallel processing.
    The code structure is inspired by Itamar Twersky's STAR_load_and_alignment_all process
    in the star.nf file from Erez Lebanon's Lab in 2023.
*/
process SECOND_STAR_MAP{
    maxForks 1
    tag "2nd MAP: ${base_comb}"
    publishDir "${params.second_map_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        path(trans_fastqs)
        path(trans_genome)
        path(trans_index_dir)

    output:
        path('*_Aligned.out*')
        // stdout

    script:
    // number of samples- for PE: number of fastq files / 2, for SE: number of files
    num_of_samples = params.pair_end ? trans_fastqs.size()/2 : trans_fastqs.size()
    // base combination (MM): first part o fthe index dir name 
    // ( for example: A2C_transformed_hg38.fa_index -> A2C)
    base_comb = trans_index_dir.name.toString().tokenize('_').get(0)
    // number of files to be mapped each time- PE: 2, SE:1 
    step = params.pair_end ? 2 : 1

    """
    # LOAD GENOME:
    echo 'loading genome'
    ${params.STAR_command} --genomeDir ${trans_index_dir} --genomeLoad LoadAndExit
    echo 'finished loading'

    # MODIFY FASQS FILES INTO A BASH LIST:  
        # transform fastq file's names into string
        fastqs_string="${trans_fastqs}"
        # Split the string into a list
        IFS=' ' read -ra samples_list <<< "\$fastqs_string"

    for element in "\${samples_list[@]}"; do
    echo "\$element"
    done

    # ALLIGNMENT:
    echo 'Starting allignmet'

    # spilitting the samples into chunks of STAR_MAX_PARALLEL
    for ((i=0; i<${num_of_samples}; i+=${params.STAR_MAX_PARALLEL})); do
        # iterate over each chunk making sure it doesn't exceed number of samples
        for ((j=i; j<i+${params.STAR_MAX_PARALLEL} && j<${num_of_samples}; j++)); do
        # extract samples and sample_id
            mate1=\${samples_list[j*${step}]}
            # SE: set mate2 to an empty string
            if [ ${params.pair_end} == 0 ]; then
                mate2=''
            # PE: set mate2 to the next file
            else
                mate2=\${samples_list[j*${step}+1]}
            fi
            sample_id=\$(echo "\${mate1}" | cut -d'_' -f2)

            echo "mate1: \${mate1}"
            echo "mate2: \${mate2}"
            echo "sample_id: \${sample_id}"
        
            # mapping each sample in the background
            ${params.STAR_command} --readFilesCommand ${params.read_files_command} --readFilesIn \${mate1} \${mate2} --genomeDir ${trans_index_dir} --outSAMattributes ${params.SAM_attr} --outSAMtype ${params.outSAMtype} --alignSJoverhangMin ${params.min_SJ_overhang} --alignIntronMax ${params.max_intron_size} --alignMatesGapMax ${params.max_mates_gap} --outFilterMismatchNoverLmax ${params.max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${params.max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${params.norm_num_of_matches} --outFilterMultimapNmax ${params.max_num_of_allignment} --genomeLoad ${params.second_map_genome_load_set} --runThreadN ${params.num_of_threads} --runDirPerm ${params.output_files_permissions} --outFileNamePrefix "./\${sample_id}_${base_comb}_" &> run_\${sample_id} &

        done

    #   Wait for the completion of each chunk's alignment to enable parallel processing,
    #   ensuring that only a maximum of STAR_MAX_PARALLEL processes are run concurrently.
        wait
    done

    ${params.STAR_command} --genomeDir ${trans_index_dir} --genomeLoad Remove
    """
}

/*
    This process utilizes an external Python script to retrieve the original sequences of the mapped reads.
    Algorithm: 
        1. Extract read IDs and original sequences from each mate's original FASTQ files (in the original_reads_dir).
        2. Iterate over each transformed mapped read in the SAM file:
            a) Check if the SAM read is mate1/mate2.
            b) Find a matching original FASTQ read from step 1.
            c) Get the original sequence and write the SAM read with the original sequence to the output file.
*/
process RETRANSFORM {
        tag "re-transform: ${sample_id} : ${base_comb}"
        publishDir "${params.retransform_output_dir}/${base_comb}", pattern: '*', mode: 'copy'

    input:
        tuple val(sample_id), path(original_reads_dir), val(base_comb), path(bam_file)
        path(python_script)

    output:
        path('*')


    script:
        outout_path = "${sample_id}_re-transformed.sam"
        """
        ${params.python_command} ${python_script} ${bam_file} ${outout_path} ${original_reads_dir}
        
        """

    
}

workflow {

    //  Get file pairs channel of Samples' both mates
    Channel
        .fromFilePairs(params.PE_reads, checkIfExists: true)
        .set {read_pairs_ch} 
    FASTP(read_pairs_ch)

    // map the quality checked fastq into the basic index dir
    unmapped_reads_ch = FIRST_STAR_MAP(FASTP.out[0], params.genome_index_dir)

    // all the bases combinations (MM) the reads to be transformed accordinly
    Channel
        .of(['A','C'],['A','G'],['A','T'],['C','A'],['C','G'],['C','T'],['G','A'],['G','C'],['G','T'],['T','A'],['T','C'],['T','G'])
        .set {bases_combination_ch}

    // Transform all of the unmapped reads recieved from FIRST_STAR_MAP (12 different bases combination)
    // group by the base combination and collect (wait) for all
    trans_reads_ch = TRANSFORM_READS(unmapped_reads_ch, bases_combination_ch)
                    // .collect(flat: false)
                    .map { it ->
                            def prefix = it[0].name.toString().tokenize('_').get(0)
                            return tuple(prefix, [it[0], it[1]])}
                    .groupTuple(by:0)
    // get all tramsformed genome and extract prefix - ref_base2alt_base, and map it to the file
    Channel
        .fromPath(params.transformed_genomes)
        .map {file -> tuple(file.name.toString().tokenize('_').get(0), file)}
        .set {transformed_genomes_ch}
    
    // transformed_genomes_ch.view()
    // get all tramsformed indexes' files
    // extract prefix - ref_base2alt_base, and map it to the files
    // make touples of [base_combination, [index_files]]
    Channel
        .fromPath(params.transformed_indexes, type: 'dir')
        .map {file -> tuple(file.name.toString().tokenize('_').get(0), file)}
        .groupTuple(by: 0)
        .set {transformed_index_ch}
    
    // concat all of the transformed files together: reads, genome, index' files
    // group by the base combination
    // for exmp: [A2C,[[A2C_read1, A2C_read2..],A2C_transformed_genome,[A2C_index_file_1, A2C_index_file_2...]]])
    // get only the files (by using map)
    tranformed_files_ch = trans_reads_ch.concat(transformed_genomes_ch)
                                        .concat(transformed_index_ch)
                                        .groupTuple(by:0)
                                        .map {tuple -> tuple[1]}

    // seperate the reads, genomes and indexes from the bases cmobination transformed files channel
    reads_channel_mount = tranformed_files_ch.map { tuple -> tuple[0].flatten() }
    genome_channel_mount = tranformed_files_ch.map { tuple -> tuple[1] }
    index_channel_mount = tranformed_files_ch.map { tuple -> tuple[2] }

    // map the transformed reads to the fitted transformed genome's index
    // collect anf faltten to get each file seperatelu and extract:
    //      sample_id (get(0)), 
    //      base_comb (get(1)).
    //      sam file
    mapped_transformed_ch = SECOND_STAR_MAP(reads_channel_mount, genome_channel_mount, index_channel_mount)
                            .collect()
                            .flatten()
                            .map { file -> tuple(
                                                file.name.toString().tokenize('_').get(0),
                                                file.name.toString().tokenize('_').get(1),
                                                file)}
    // get the dirs of the original reads (output of the first map process) and extract:
    //      sample id
    //      file (directory)
    Channel
        .fromPath(params.original_reads, type:'dir')
        .map {file -> tuple(file.name.toString(), file)}
        .set {originial_reads_ch}

    // combine the mapped transformed sam files with the original fastqs using the sample id as key
    // and retransform the sequences of the mapped bam to the originak sequences
    RETRANSFORM(originial_reads_ch.combine(mapped_transformed_ch, by:0), params.retransform_python_script)
    RETRANSFORM.out.view()
 
}

