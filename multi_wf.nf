process GILI {
    output:
        stdout
    script:
    """
    echo gili 
    """
}

process SILI {
    output:
        stdout
    script:
    """
    echo sili 
    """
}

workflow gili {
    gili_ch = GILI()
    gili_ch.view()

}

workflow {
    sili_ch = SILI()
    sili_ch.view()
}