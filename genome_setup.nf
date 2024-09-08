"""
author: Gili Wolf 29-11-2023
---------------------------------------
HE-GENOME SETUP:
this scripts:
    1.   changes each refrence base in the original genome
          into a new alternative base, in order to create transformed genome.
    2.   create genome indexes by using STAR --runMode genomeGenerate.
As part of the HyperEditing tool, this is the pre-setup to further analysis in the 2nd mapping stage. 

usage: 
nextflow genome_setup.nf -c genome_setup.nf.config --genome_fasta <PATH_TO_GENOME> --genome_setup_outdir <PATH_TO_OUTOUT_DIR>

"""

bases=['A','G','C','T']
def helpMessage() {
    log.info"""
            HYPER-EDITING GENOME-SETUP
            ==========================
            transform and index ${params.genome_fasta} 12 times as follows:
            a->g
            a->c
            a->t

            g->a
            g->c
            g->t

            c->a
            c->g
            c->t

            t->a
            t->g
            t->c
    """.stripIndent()
}
// transform refrence bases to alternative bases in the fasta files into a new file
// in generic_transform/transformed_genome
process TRANSFORM {
    tag "${ref_base} >> ${alt_base}"
    maxForks 5
    publishDir "${params.transform_genome_output_dir}", pattern: '*.{fa,fasta}', mode: 'copy'

    input:
        each ref_base
        each alt_base
        path ref_genome
    
    output:
        stdout
        path('*.{fa,fasta}')


    // condition: don't transform if refrence base and alterantive base are the same
    when:
        ref_base != alt_base

    shell:
        '''
        
        # output files:
        transformed_genome_path="!{ref_base}2!{alt_base}_transformed_!{ref_genome}"
        index_genome_path="!{ref_base}2!{alt_base}_transformed_!{ref_genome}_index"


        # transform each letter in the genome from ref_bas > alt_base
        transform_genome() {
            #create lowercase ref and alt bases
            ref_base_lower=$(echo "!{ref_base}" | tr '[:upper:]' '[:lower:]')
            alt_base_lower=$(echo "!{alt_base}" | tr '[:upper:]' '[:lower:]')
            #create uppercase ref and alt bases
            ref_base_upper=$(echo "!{ref_base}" | tr '[:lower:]' '[:upper:]')
            alt_base_upper=$(echo "!{alt_base}" | tr '[:lower:]' '[:upper:]')


            # transform all lowercase ref bases into lowercase alt bases
            # transform all uppercase ref bases into uppercase alt bases
            awk -v ref_base_lower="${ref_base_lower}" -v alt_base_lower="${alt_base_lower}" \
            -v ref_base_upper="${ref_base_upper}" -v alt_base_upper="${alt_base_upper}" \
            '/^>/ {print; next} {gsub(ref_base_lower, alt_base_lower); gsub(ref_base_upper, alt_base_upper); print}' \
            "!{ref_genome}" > "${transformed_genome_path}"
        
        }

        #transformation
            echo -e "---------------------\nstarted tranforming ${transformed_genome_path}\n"
            transform_genome
            echo -e "********************\nended transforming ${transformed_genome_path}\n---------------------\n"

        '''
}

// indexing all of the transformed genome using STAR (from docker)
process INDEX() {
    tag "${transformed_fasta}"
    maxForks 1
    // container 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0'
    publishDir "${params.index_output_dir}", pattern: '*_index', mode: 'copy' //specific index dir
    publishDir "${params.index_output_dir}/${transformed_fasta}_index", pattern: '*[!_index]', mode: 'copy' //index files
    

    input:
        path transformed_fasta

    output:
        stdout
        path('*')

    shell:
    '''
        # indexing dir
        index_genome_dir="!{transformed_fasta}_index"

        # index the transformed genome using STAR genomeGenerate
        star_genome_indexing() {
            local star_program='STAR'
            local num_of_threads=20
            local index_command="${star_program} --runThreadN ${num_of_threads} --runMode genomeGenerate --genomeDir ${index_genome_dir} --genomeFastaFiles !{transformed_fasta} --runDirPerm All_RWX"

            mkdir ${index_genome_dir}

            #evaluate the index command
            ${index_command}
            
        }

            #indexing
            echo -e "---------------------\nstarted indexing ${index_genome_dir}\n"
            star_genome_indexing
            echo -e "********************\nended indexing ${index_genome_dir}\n---------------------\n"

            '''

}

workflow GENOME_SETUP {
    take: 
        genome_ch
    main:
        // 12 transformation of the genome
        TRANSFORM(bases,bases,genome_ch)
    
        // indexing of the 12 transformed genomes
        INDEX(TRANSFORM.out[1])
    emit:
        INDEX.out[1]
}


workflow {
    if(params.help){
        helpMessage()
        System.exit()
    }
    // genome fasta channel
    Channel
        .fromPath(params.genome_fasta, checkIfExists: true)
        .set {genome_ch}
    
    GENOME_SETUP(genome_ch)
}   




    