// file settings
docker {
    enabled = true
    docker.runOptions='-u $(id -u):$(id -g)'
}

executor {
    $local {
        queueSize = 10
  }
}

params {
    // Input paths:
    SE_HE_reads="/private10/Projects/Gili/HE_workdir/first_part/SE_test/re-transform/**/*.sam"
    pair_end=0
    fasta_path="/private10/Projects/Gili/HE_workdir/genome_setup/hg38.fa"
    detect_python_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/detect_clusters.py"
    filter_python_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/filter_clusters.py"
    file_seperator="_"
    python_command ="python"

    // Output dir:
    outdir = "detection/first_try"
    detect_output_dir = "${params.outdir}/detected_clusters"
    filter_output_dir = "${params.outdir}/filtered_clusters"

}

// container settings:
process {
    withName: DETECT {
        container = 'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0'
    }
    withName: FILTER {
        container = 'quay.io/jupyter/scipy-notebook'
    }
}

