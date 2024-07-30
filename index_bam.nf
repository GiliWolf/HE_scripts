process INDEX_BAM {
        maxForks 1
        tag "index: ${sample_id}"
        publishDir "/private10/Projects/Gili/HE_workdir/detection/INDEX_TRY", pattern: '*', mode: 'copy'

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path(bam_file), path("*.bai"), env(COUNT)

    script:
        sorted_bam_path = "${sample_id}_sorted.bam"

        """
        samtools sort -o ${sorted_bam_path} ${bam_file}
        samtools index -@ 16 ${sorted_bam_path}
        COUNT=\$(samtools view -c ${sorted_bam_path})

        """ 
}

workflow {

    // GET SAMPLES:
    Channel
        .fromPath("/private10/Projects/Gili/HE_workdir/first_part/GTEX/ArteryAorta_PE/second_map/A2C/GTEX-RUSQ_ArteryAorta_GTEX-RUSQ-0326-SM-47JWS_A2C_Aligned.out.bam", checkIfExists: true)
        .map {file -> tuple (file.baseName, file)}
        .set {samples_ch} 

    index_ch = INDEX_BAM(samples_ch)
    index_ch.view()
}