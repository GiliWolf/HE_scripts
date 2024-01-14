process TEST {
    container = 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0'
    input:
        path samples

    output:
        stdout
        // path(*)
      script:
    for sample in ${samples} {
        """
        STAR --help
        """
    }
  
}

workflow{
    Channel
        .fromPath("/private10/Projects/Gili/HE_workdir/first_part/SRA2/*.fastq")
        .set {samples_ch}
    
    TEST(samples_ch.collect())
    TEST.out.view()
}
