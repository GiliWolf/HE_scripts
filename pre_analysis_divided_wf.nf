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
params.help=false

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

def helpMessage() {
    log.info '''
        usage: nextflow -c pre_anylis.nf.config [options...] run pre_anlysis.nf
        * for running in baskground:
          nohup nextflow -c pre_anylis.nf.config [options...] -bg run pre_anlysis.nf &
        * for help message: nextflow -c run pre_anlysis.nf --help

        options:
        ---------------
                --help              displayes help message             
        ---------------
        input file parameters 
                -pair_end			pair end flag (0:SE, 1:PE)
                -reads_dir= 			samples' directory
                -reads_suffix=			suffixs of the samples' files (default-".fastq")
                -file_seperator=		seperator of the file's name companats (default- '_')
                -mate_seperator			suffix of both mates (default- "_")
                -mate1_suff			suffix of mate1 (default- "1")
                -mate2_suff			suffix of mate2 (default- "2")
        ---------------
        output dirs parameters
                -outdir				path for the output directory (relative path preferred)
                -fastp_output_dir		path for the output of the fastp process (default- "${params.outdir}/fastp")
                -first_map_output_dir		path for the output of the first map process
                                                (default-"${params.outdir}/first_map")
                -transform_output_dir		path for the output of the transform reads process 
                                                (default-"${params.outdir}/transformed_unmapped")
                -second_map_output_dir		path for the output of the second map process
                                                (default- "${params.outdir}/second_map")
                -retransform_output_dir		path for the output of the second map process
                                                (default- "${params.outdir}/re-transform")
        ---------------
        Fastp parameters
                -N_bases_num			maximum number of N bases in a read (default- 5)
                -avg_quality			minimum average quality score of a read (default- 30)
                -low_quality_per 		minimum percentage of bases allowed to be unqualified (default- 20)
                -low_quality_num		minimum quality value that a base is qualified, Phred score (default- 25)
        ---------------
        first map parameters
                -fastp_output_suffix		suffix of the unmapped files (default - ".processed.fastq")
                -read_files_command		command line to execute for each of the input file (default - "cat")
                -genome_index_dir		path to the genome's index directory
                -SAM_attr			desired SAM attributes (default - "All")
                -outSAMtype			type of the output mapped files (default - "BAM Unsorted")
                -min_SJ_overhang		minimum overhang (i.e. block size) for spliced alignments (default -8)
                -max_intron_size		maximum intron length (dafault value for STAR and HE- 1000000)
                -max_mates_gap			maximum genomic distance between mates (default -600000)
                -max_mismatches_ratio_to_ref	max ratio of mismatches to *mapped* length (dafault value for STAR and HE- 3)
                -max_mismatche_ratio_to_read	max ratio of mismatches to *read* length (dafault value for STAR and HE- 1)
                -norm_num_of_matches		min number of matched bases normalized to the read length
                                                (sum of mates’ lengths for PE reads) (dafault value for STAR and HE- 0.66)
                -max_num_of_allignment		max number of multiple alignments allowed for a read, if exceeded->unmapped (default - 1)
                -genome_load_set		genome shared memory is not used (dafault value for STAR and HE-"NoSharedMemory")
                -num_of_threads			number of threads for each mapping process (dafault - 5)
                -unmapped_out_files		unmapped reads will be output into separate file(s)Unmapped.out.mate1[2],
                                            formatted the same way as -input read files (default - Fastx")
                -output_files_permissions	file permissions (dafault - "All_RWX")
        ---------------
        second map parameters
        general:
                -transformed_indexes		path to the directory of the transformed genome indexes from HE part 1
                                                (for example - "/generic_transform/hg38transform/genome_indexing/*")
        STAR parameters:
                -STAR_MAX_PARALLEL		number of files to be mapped parallely (dafault - 6)
                -second_map_genome_load_set	controls how the genome is loaded into memory (dafault -"LoadAndKeep")
        ---------------
        retransform parameters
                -original_reads			path to the original fastq files (default - extracts from first_map_output_dir)
                -filter_sam_files		STAR output files to be output (default - '*Aligned.out*')
    '''
    .stripIndent()
}

/*
    Preprocess the reads using FASTP to filter out low-quality reads.
*/
process FASTP{
    tag "preprocessing: ${sample_id}"
    publishDir "${params.fastp_output_dir}", pattern: '*.{fastq,json}', mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        path('*.fastq')
        path('*.json')


    script:

        if (params.pair_end == 0)
            """
            # SE (run on all files ending with the reads suffixs)
                        ${params.fastp_command} \
                        -n ${params.N_bases_num} \
                        -e ${params.avg_quality} \
                        -u ${params.low_quality_per} \
                        -q ${params.low_quality_num} \
                        --low_complexity_filter \
                        --complexity_threshold ${params.complexity_threshold} \
                        -j "${sample_id}.fastp.json" \
                        --dont_eval_duplication \
                        --in1 "${reads}" \
                        -o "${sample_id}${params.fastp_output_suffix}"
            """
        
        else
            """
                #PE (run on all 1st mate files, and change the mate suffixs for 2nd mate)
                ${params.fastp_command} -n ${params.N_bases_num} -e ${params.avg_quality} -u ${params.low_quality_per} -q ${params.low_quality_num} --low_complexity_filter --complexity_threshold ${params.complexity_threshold} -j "${sample_id}.fastp.json" --dont_eval_duplication --in1 ${reads[0]} -o "${sample_id}${params.mate_seperator}${params.mate1_suff}${params.fastp_output_suffix}" --in2 ${reads[1]} -O "${sample_id}${params.mate_seperator}${params.mate2_suff}${params.fastp_output_suffix}"
            """

}

/*
    Map all the reads using STAR with the original genome's index to identify the unmapped reads for further processing.
*/
process FIRST_STAR_MAP{
    maxForks 1
    tag "1st MAP: ${genome_index}"
    publishDir "${params.first_map_output_dir}", pattern: "*", mode: 'copy'

    input:
        path(fastq)
        path(genome_index)

    output:
        path('*Unmapped*'), emit: fastqs
        path('*bam')

    script:
    // number of samples- for PE: (number of fastq files / 2), for SE: number of files
    num_of_samples = params.pair_end ? fastq.collect().size()/2 : fastq.collect().size()
    // number of files to be mapped each time- PE: 2, SE:1 
    step = params.pair_end ? 2 : 1


    script:
    """
    # LOAD GENOME:
    echo 'loading genome'
    ${params.STAR_command} --genomeDir ${genome_index} --genomeLoad LoadAndExit
    echo 'finished loading'

    # MODIFY FASTQS FILES INTO A BASH LIST:  
        # transform fastq file's names into string
        fastqs_string="${fastq}"
        # Split the string into a list
        IFS=' ' read -ra samples_list <<< "\$fastqs_string"

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

            # extract sample ID - same logic as to the groovy function getSampleID (above the main workflow)
            # only it is written in bash
            if [ "${params.pair_end}" -eq 0 ]; then
                sample_id=\$(echo "\${mate1}" | awk -v fs="${params.file_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ss}{split(\$1, parts, fs); for(i=1; i<length(parts); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)]}')
            else
                if [[ "${params.mate_seperator}" == "${params.file_seperator}" ]]; then
                    sample_id=\$(echo "\${mate1}" | awk -v fs="${params.file_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ss}{split(\$1, parts, fs); for(i=1; i<(length(parts)-1); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)-1]}')
                else
                    sample_id=\$(echo "\${mate1}" | awk -v ms="${params.mate_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ms ss}{print \$1}')
                fi

            fi

        
            # mapping each sample in the background
            ${params.STAR_command} --readFilesCommand ${params.read_files_command} --readFilesIn \${mate1} \${mate2} --genomeDir ${genome_index} --outSAMattributes ${params.SAM_attr} --outSAMtype ${params.outSAMtype} --alignSJoverhangMin ${params.min_SJ_overhang} --alignIntronMax ${params.max_intron_size} --alignMatesGapMax ${params.max_mates_gap} --outFilterMismatchNoverLmax ${params.max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${params.max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${params.norm_num_of_matches} --outFilterMultimapNmax ${params.max_num_of_allignment_first_map} --genomeLoad ${params.genome_load_set} --runThreadN ${params.num_of_threads} --outReadsUnmapped ${params.unmapped_out_files} --runDirPerm ${params.output_files_permissions} --outFileNamePrefix "./\${sample_id}." &> run_\${sample_id} &

        done

    #   Wait for the completion of each chunk's alignment to enable parallel processing,
    #   ensuring that only a maximum of STAR_MAX_PARALLEL processes are run concurrently.
        wait
    done
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

        if (params.pair_end == 0)   //SE
            """
                # create lowercase ref and alt bases
                    ref_base_lower=\$(echo "${ref_base}" | tr '[:upper:]' '[:lower:]')
                    alt_base_lower=\$(echo "${alt_base}" | tr '[:upper:]' '[:lower:]')

                # create uppercase ref and alt bases
                    ref_base_upper=\$(echo "${ref_base}" | tr '[:lower:]' '[:upper:]')
                    alt_base_upper=\$(echo "${alt_base}" | tr '[:lower:]' '[:upper:]')
                
                # create output path
                    transformed_output_path="${ref_base}2${alt_base}${params.file_seperator}${sample_id}.fastq"
            
                # generic awk script for each mate
                    awk_script='
                        NR%4==2 {
                            gsub(ref_base_lower, alt_base_lower);
                            gsub(ref_base_upper, alt_base_upper)
                        } 1
                    '

                # transform
                    awk -v ref_base_lower="\${ref_base_lower}" -v alt_base_lower="\${alt_base_lower}" \
                    -v ref_base_upper="\${ref_base_upper}" -v alt_base_upper="\${alt_base_upper}" \
                    "\${awk_script}" "${unmapped_fastq}" > "\${transformed_output_path}"
            """

        else                        //PE
            """
                # create lowercase ref and alt bases
                    ref_base_lower=\$(echo "${ref_base}" | tr '[:upper:]' '[:lower:]')
                    alt_base_lower=\$(echo "${alt_base}" | tr '[:upper:]' '[:lower:]')

                # create uppercase ref and alt bases
                    ref_base_upper=\$(echo "${ref_base}" | tr '[:lower:]' '[:upper:]')
                    alt_base_upper=\$(echo "${alt_base}" | tr '[:lower:]' '[:upper:]')
                
                # create output path for each mate
                    transformed_output_path_mate1="${ref_base}2${alt_base}${params.file_seperator}${sample_id}_1.fastq"
                    transformed_output_path_mate2="${ref_base}2${alt_base}${params.file_seperator}${sample_id}_2.fastq"
            
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
        path(trans_index_dir)

    output:
        path('*_Aligned.out*')
        // stdout

    script:
    // number of samples- for PE: (number of fastq files / 2), for SE: number of files
    num_of_samples = params.pair_end ? trans_fastqs.collect().size()/2 : trans_fastqs.collect().size()
    // base combination (MM): first part o fthe index dir name 
    // ( for example: A2C_transformed_hg38.fa_index -> A2C)
    base_comb = trans_index_dir.name.toString().tokenize(params.file_seperator).get(0)
    // number of files to be mapped each time- PE: 2, SE:1 
    step = params.pair_end ? 2 : 1

    """
    # LOAD GENOME:
    echo 'loading genome'
    ${params.STAR_command} --genomeDir ${trans_index_dir} --genomeLoad LoadAndExit
    echo 'finished loading'

    # MODIFY FASTQS FILES INTO A BASH LIST:  
        # transform fastq file's names into string
        fastqs_string="${trans_fastqs}"
        # Split the string into a list
        IFS=' ' read -ra samples_list <<< "\$fastqs_string"


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

        # extract sample ID - *similar* to the logic of the groovy function getSampleID (above main workflow)
        # *execpt* of removing the BASE_COMB which located in the beginning of the file.
        # for example - "A2C_ADAR_GMCSF_AdarWT_MDA5KO_79_372_SKO-2.fq" -> "ADAR_GMCSF_AdarWT_MDA5KO_79_372_SKO"
        if [ "${params.pair_end}" -eq 0 ]; then
            sample_id=\$(echo "\${mate1}" | awk -v fs="${params.file_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ss}{split(\$1, parts, fs); for(i=2; i<length(parts); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)]}')
        else
            if [[ "${params.mate_seperator}" == "${params.file_seperator}" ]]; then
                sample_id=\$(echo "\${mate1}" | awk -v fs="${params.file_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ss}{split(\$1, parts, fs); for(i=2; i<(length(parts)-1); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)-1]}')
            else
                sample_id=\$(echo "\${mate1}" | awk -v fs="${params.file_seperator}" -v ms="${params.mate_seperator}" -v ss="${params.suffix_seperator}" 'BEGIN{FS=ms ss}{split(\$1, parts, fs); for(i=2; i<length(parts); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)]}')
            fi
        fi


            # mapping each sample in the background
            ${params.STAR_command} --readFilesCommand ${params.read_files_command} --readFilesIn \${mate1} \${mate2} --genomeDir ${trans_index_dir} --outSAMattributes ${params.SAM_attr} --outSAMtype ${params.outSAMtype} --alignSJoverhangMin ${params.min_SJ_overhang} --alignIntronMax ${params.max_intron_size} --alignMatesGapMax ${params.max_mates_gap} --outFilterMismatchNoverLmax ${params.max_mismatches_ratio_to_ref} --outFilterMismatchNoverReadLmax ${params.max_mismatche_ratio_to_read} --outFilterMatchNminOverLread  ${params.norm_num_of_matches} --outFilterMultimapNmax ${params.max_num_of_allignment_second_map} --genomeLoad ${params.second_map_genome_load_set} --runThreadN ${params.num_of_threads} --runDirPerm ${params.output_files_permissions} --outFileNamePrefix "./\${sample_id}${params.file_seperator}${base_comb}${params.file_seperator}" &> run_\${sample_id} &

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
        tuple val(sample_id), path(original_reads), val(base_comb), path(bam_file)
        path(python_script)

    output:
        path('*'), optional: true

    // python re-transform.py -s <sam_file> -o <output_path> -i <original_fastq_path> [-I <original_fastq_path_2>] <pair_end[0/1]>
    script:
        outout_path = "${base_comb}${params.file_seperator}${sample_id}${params.file_seperator}re-transformed.sam"
        if (params.pair_end == 0)   //SE
            """
            ${params.python_command} ${python_script} -s ${bam_file} -o ${outout_path} -i ${original_reads} -pe ${params.pair_end}
            
            """
        else
            """
            ${params.python_command} ${python_script} -s ${bam_file} -o ${outout_path} -i ${original_reads[0]} -I ${original_reads[1]} -pe ${params.pair_end}
            """

}


workflow pair_end {
    take:
        samples_ch
        bases_combination_ch
        transformed_index_ch

    main:
        FASTP(samples_ch)
        fastp_ch = FASTP.out[0].collect()

        // if ( params.flag ) {
        // bar()
        // omega_ch = bar.out
        // }
        // else {
        //     foo()
        //     omega_ch = foo.out
        // }
        //
        FIRST_STAR_MAP(fastp_ch, params.genome_index_dir)
        FIRST_STAR_MAP.out.fastqs
                        .flatten()
                        .map { file ->
                                def file_sample_id  = file.name.toString().tokenize(params.suffix_seperator).get(0)
                                return tuple(file_sample_id,file)
                        }
                        .groupTuple(by:0)
                        .set {unmapped_reads_ch}
        //
        trans_reads_ch = TRANSFORM_READS(unmapped_reads_ch, bases_combination_ch)
                        .map { file ->
                                def prefix = file[0].name.toString().tokenize(params.file_seperator).get(0)
                                return tuple(prefix, [file[0], file[1]])}
                        .groupTuple(by:0)
        tranformed_files_ch = trans_reads_ch
                                    .concat(transformed_index_ch)
                                    .groupTuple(by:0)
                                    .map {tuple -> tuple[1]}
        reads_channel_mount = tranformed_files_ch.map { tuple -> tuple[0].flatten() }
        index_channel_mount = tranformed_files_ch.map { tuple -> tuple[1] }
        //
        mapped_transformed_ch = SECOND_STAR_MAP(reads_channel_mount, index_channel_mount)
                                .collect()
                                .flatten()
                                .map { file -> tuple(
                                file.name.toString().split(/_[A-Z]2[A-Z]_/)[0],
                                file.name.toString().tokenize(params.file_seperator)[-2],
                                file)} 
        
        // Channel
        //     .fromPath("/private10/Projects/Gili/HE_workdir/first_part/GTEX/MuscleSkeletal_PE/second_map/**/*_Aligned.out*")
        //     .collect()
        //     .flatten()
        //     .map { file -> tuple(
        //                         file.name.toString().split(/_[A-Z]2[A-Z]_/)[0],
        //                         file.name.toString().tokenize(params.file_seperator)[-2],
        //                         file)} 
        //     .set{mapped_transformed_ch}
        // Channel
        //     .fromFilePairs(params.PE_original_reads)
        //     .set {originial_reads_ch}
        //
        // unmapped_reads_ch.view()
        // mapped_transformed_ch.view()
        files_to_retransform_ch = unmapped_reads_ch.combine(mapped_transformed_ch, by:0)
        // files_to_retransform_ch.view()
        retransform_ch = RETRANSFORM(files_to_retransform_ch, params.retransform_python_script)

    // emit:
}

workflow single_end {
    take:
        samples_ch
        bases_combination_ch
        transformed_index_ch

    main:
        FASTP(samples_ch)
        fastp_ch = FASTP.out[0].collect()
        //
        FIRST_STAR_MAP(fastp_ch, params.genome_index_dir)
        FIRST_STAR_MAP.out.fastqs
                        .flatten()
                        .map { file ->
                                def file_sample_id  = file.name.toString().tokenize(params.suffix_seperator).get(0)
                                return tuple(file_sample_id,file)
                        }
                        .groupTuple(by:0)
                        .set {unmapped_reads_ch}
        //
        trans_reads_ch = TRANSFORM_READS(unmapped_reads_ch, bases_combination_ch)
                        .map { file ->
                                def prefix = file.name.toString().tokenize(params.file_seperator).get(0)
                                return tuple(prefix, file)}
                        .groupTuple(by:0)
        tranformed_files_ch = trans_reads_ch
                                    .concat(transformed_index_ch)
                                    .groupTuple(by:0)
                                    .map {tuple -> tuple[1]}
        reads_channel_mount = tranformed_files_ch.map { tuple -> tuple[0].flatten() }
        index_channel_mount = tranformed_files_ch.map { tuple -> tuple[1] }
        //
        mapped_transformed_ch = SECOND_STAR_MAP(reads_channel_mount, index_channel_mount)
                                .collect()
                                .flatten()
                                .map { file -> tuple(
                                file.name.toString().split(/_[A-Z]2[A-Z]_/)[0],
                                file.name.toString().tokenize(params.file_seperator)[-2],
                                file)} 
        //
        Channel
            .fromPath(params.original_reads)
            .map {file -> tuple(file.name.toString().tokenize(params.suffix_seperator).get(0), file)}
            .set {originial_reads_ch}
        //
        files_to_retransform_ch = originial_reads_ch.combine(mapped_transformed_ch, by:0)
        files_to_retransform_ch.view()
        retransform_ch = RETRANSFORM(files_to_retransform_ch, params.retransform_python_script)


}


workflow {

    if(params.help){
        helpMessage()
        System.exit(1)
    }

    /*
    *   GET SAMPLES:
    */ 
    if (params.pair_end == 0)                         //SE
        Channel
            .fromPath(params.SE_reads, checkIfExists: true)
            .map {file -> tuple (file.baseName, file)}
            .set {samples_ch}
    else if (params.pair_end == 1)                    //PE
        Channel
            .fromFilePairs(params.PE_reads, checkIfExists: true)
            .set {samples_ch} 
    else             // raise error if pair_end flag != 0/1
        error "----------------\n error: pair_end flag must be 0/1\n----------------"
    
    //
    Channel
        .of(['A','C'],['A','G'],['A','T'],['C','A'],['C','G'],['C','T'],['G','A'],['G','C'],['G','T'],['T','A'],['T','C'],['T','G'])
        .set {bases_combination_ch}
    //
    Channel
        .fromPath(params.transformed_indexes, type: 'dir')
        .map {file -> tuple(file.name.toString().tokenize(params.file_seperator).get(0), file)}
        .groupTuple(by: 0)
        .set {transformed_index_ch}


    if (params.pair_end == 0)  
        single_end(samples_ch, bases_combination_ch, transformed_index_ch)
    else
        pair_end(samples_ch, bases_combination_ch, transformed_index_ch)

}


// output {
//     directory params.top_outdir
// }


