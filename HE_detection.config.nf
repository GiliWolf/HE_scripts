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
    grid_serch_python_script = "/private10/Projects/Gili/HE_workdir/HE_scripts/HE_grid_search.py"
    file_seperator="_"
    python_command ="python"

    // Output dir:
    outdir = "detection/with_grid"
    detect_output_dir = "${params.outdir}/detected_clusters"
    filter_output_dir = "${params.outdir}/filtered_clusters"

    // detect script argumnets - 
    columns_select = 'all'
    // ranges of grid search
    // editing sites - 
    es_start = 0
    es_end = 1
    es_step = 1

    // editing fracture
    ef_start = 0.5
    ef_end = 0.6
    ef_step = 0.1

    
    // phred score
    ps_start = 30
    ps_end = 40
    ps_step = 10

    
    // num of es to read length ratio
    es2l_start = 0.03
    es2l_end = 0.04
    es2l_step = 0.01

    // cluster length to read length ratio
    cl2l_start = 0.01
    cl2l_end = 0.02
    cl2l_step = 0.01




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

