process CPAT_PREDICT {
    tag "cpat"
    label 'process_low'
    // Decision-record D3: CPAT 3.0.5 replaces the python2 CNCI (retired with
    // the evidence recorded in the SDR). Species models are params: the
    // human files ship with CPAT; a rice model is user-provided.
    container 'quay.io/biocontainers/cpat:3.0.5--py310h20b60a1_4'

    input:
    tuple val(meta), path(fasta)
    path model
    path hexamer

    output:
    tuple val(meta), path("cpat.out.tsv"), emit: tsv

    script:
    def args = task.ext.args ?: ''
    """
    cpat.py \\
        -c ${model} \\
        -d ${hexamer} \\
        -f ${fasta} \\
        -o cpat.out.tsv \\
        -g ${meta.id} \\
        ${args}
    """

    stub:
    """
    touch cpat.out.tsv
    """
}
