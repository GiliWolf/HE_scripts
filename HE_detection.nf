process DETECT {
        tag "detection: ${sample_id}"
        publishDir "", pattern: '*', mode: 'copy'

    input:


    output:
        path('*')

    script:

        """
        
        """  
}

process FILTER {
        tag "filter: ${sample_id}"
        publishDir "", pattern: '*', mode: 'copy'

    input:


    output:
        path('*')

    script:

        """
        
        """  
}

workflow {  
    DETECT()

    FILTER()

}