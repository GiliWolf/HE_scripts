
include { GENOME_SETUP } from '/private10/Projects/Gili/HE_workdir/HE_scripts/genome_setup.nf' 
include { PRE_ANALYSIS } from '/private10/Projects/Gili/HE_workdir/HE_scripts/pre_analysis_divided_wf.nf'
include { HE_DETECTION } from '/private10/Projects/Gili/HE_workdir/HE_scripts/HE_detection.nf'
// TRANSFORM: 
// ./nextflow HE_scripts/MAIN_HE_SCRIPT.nf -c HE_scripts/MAIN_HE_SCRIPT.nf.config --run_mode {transformGenome | align | detect | align-detect} [TRANSFORM: --genome_fasta <PATH_TO_GENOME> --genome_setup_outdir MAIN_tests/genome_setup] | [ALIGN: --reads_dir <READS_FATQ_PATH> --outdir MAIN_tests/pre_analysis --genome_index_dir genome_setup/hg38_index --transform_genome_dir genome_setup/generic_transform/hg38transform --pair_end {0,1}] | [DETECT: --detect_input_dir <PATH_TO_ALIGN_OUTPUT_DIR> --detect_outdir MAIN_tests/detect --fasta_path "/private10/Projects/Gili/HE_workdir/genome_setup/hg38.fa --fasta_index_path "/private10/Projects/Gili/HE_workdir/genome_setup/hg38.fa.fai"]

def helpMessage() {
    log.info '''
        MAIN HE SCRIPT: 
        
    '''
    .stripIndent()
}

workflow transformGenome {
    main:
        Channel
            .fromPath(params.genome_fasta, checkIfExists: true)
            .set {genome_ch}
    
        GENOME_SETUP(genome_ch)
}

workflow align {
    main:
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
        
        potential_HE_ch = PRE_ANALYSIS(samples_ch)
    emit:
        potential_HE_ch
}

workflow detect {
    main: 
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
                .map {file -> tuple (file.name.toString().tokenize(params.file_seperator).get(0),file.name.toString().tokenize(params.suffix_seperator).get(0), file)}
                .groupTuple(by: [0,1], sort:true)
                .map {tuple -> tuple.flatten()}
                .set {samples_ch}
        }
        samples_ch = COUNT_RECORDS(samples_ch)
        HE_DETECTION(samples_ch)
}

workflow detect_from_align {
    take:
        bam_files_ch
    
    main:
        HE_DETECTION(bam_files_ch)

}
workflow {

    if(params.help){
        helpMessage()
        System.exit()
    }

    //GENOME SETUP: TRANSFORM, INDEX
    if (params.run_mode == "transformGenome"){
        transformGenome()
    }
    //PRE-ANALYSIS: FATSP, 1ST MAP, TRANSFORM, 2ND MAP, RETRANSFORM
    else if (params.run_mode == "align") {
        align()
    }
    // DETECTION: INDEX, DETECT, FILTER
    else if (params.run_mode == "detect"){
        detect()
    }
    // SERACH & DETECT
    else if (params.run_mode == "align-detect"){
        potential_HE_ch = align()
        detect_from_align(potential_HE_ch)
    }
    else {
        error "HE usage: --run_mode {transformGenome | align | detect | align-detect}"
    }
}