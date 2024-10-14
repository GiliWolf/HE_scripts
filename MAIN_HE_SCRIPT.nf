
include { GENOME_SETUP } from '/private10/Projects/Gili/HE_workdir/HE_scripts/genome_setup.nf' 
include { PRE_ANALYSIS } from '/private10/Projects/Gili/HE_workdir/HE_scripts/pre_analysis_divided_wf.nf'
include { HE_DETECTION } from '/private10/Projects/Gili/HE_workdir/HE_scripts/HE_detection.nf'

/*
----------------------------
Author: Gili Wolf
Date: 01-08-2024
Affiliation: Levanon Lab
----------------------------

Description:
    This script controls the workflow of the Hyper-Editing detection tool, controlling of three key steps:

    1. Genome Transformation:
        The script performs genome transformations, creating 12 different transformations for each possible base combination, and indexes the transformed genomes using STAR. This step generates a transformed genome that facilitates the detection of hyper-editing events.

    2. Pre-Analysis:
        This step is the core of the original hyper-editing pipeline. It aligns the sample reads to the original genome, then transforms unmapped reads and maps them again to the transformed genome. The original sequences of the transformed reads from the second alignment are then recovered and stored in alignment BAM files.

    3. Detection:
        This step involves analyzing the BAM files to detect editing sites and mismatches, filtering the results to identify hyper-edited reads.

Run Modes:
    The script offers multiple run modes to customize the workflow based on user needs:
    - transformGenome: Executes only the Genome Transformation step.
    - align: Executes only the Pre-Analysis step.
    - detect: Executes only the Detection step.
    - align-detect: Combines Pre-Analysis and Detection steps, running both consecutively

Usage:

    ./nextflow -c HE_scripts/MAIN_HE_SCRIPT.nf.config run HE_scripts/MAIN_HE_SCRIPT.nf --
    run_mode {transformGenome|align|detect|align-detect} [TRANSFORM: --genome_fasta <PATH_TO_GENOME> --genome_setup_outdir <OUTDIR_PATH>] | [ALIGN: --align_reads_dir <READS_FATQ_PATH> --align_outdir <OUTDIR_PATH> --genome_index_dir <GENOME_INDEX_PATH> --transform_genome_dir <TRANSFORM_GENOME_DIR_PATH> --pair_end {0,1}] | [DETECT: --detect_input_dir <PATH_TO_ALIGN_OUTPUT_DIR> --detect_outdir <OUTDIR_PATH> --genome_fasta <PATH_TO_GENOME> --genome_index_dir <GENOME_INDEX_PATH>]

*/

params.help=false

def helpMessage() {
    log.info '''
        HYPER-EDITING MAIN SCRIPT:
        ===========================
        This script controls the workflow of the Hyper-Editing detection tool.
        Avaiable run modes:
        - transformGenome: Executes only the Genome Transformation step.
        - align: Executes only the Pre-Analysis (main algorithm) step.
        - detect: Executes only the Detection step.
        - align-detect: Combines Pre-Analysis and Detection steps, running both consecutively
        ===========================
        usage: ./nextflow -c HE_scripts/MAIN_HE_SCRIPT.nf.config run HE_scripts/MAIN_HE_SCRIPT.nf --
        run_mode {transformGenome|align|detect|align-detect} [TRANSFORM: --genome_fasta <PATH_TO_GENOME> --genome_setup_outdir <OUTDIR_PATH>] | [ALIGN: --align_reads_dir <READS_FATQ_PATH> --align_outdir <OUTDIR_PATH> --genome_index_dir <GENOME_INDEX_PATH> --transform_genome_dir <TRANSFORM_GENOME_DIR_PATH> --pair_end {0,1}] | [DETECT: --detect_input_dir <PATH_TO_ALIGN_OUTPUT_DIR> --detect_outdir <OUTDIR_PATH> --genome_fasta <PATH_TO_GENOME> --genome_index_dir <GENOME_INDEX_PATH>]
        
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