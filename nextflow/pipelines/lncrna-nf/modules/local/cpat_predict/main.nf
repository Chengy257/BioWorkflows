process CPAT_PREDICT {
    tag "cpat"
    label 'process_low'
    // Decision-record D3: CPAT 3.0.5 replaces the python2 CNCI (retired with
    // the evidence recorded in the SDR). Species models are params: the
    // human files ship with CPAT; a rice model is user-provided.
    // CPAT 3.x CLI: executable 'cpat'; -g FASTA, -d model, -x hexamer, -o is
    // an output PREFIX (not a file name). It writes <prefix>.ORF_prob.best.tsv
    // (best ORF per sequence, Coding_Probability column) and, when sequences
    // have no ORF at all, <prefix>.no_ORF.txt (CPAT considers these
    // non-coding). The touch keeps the output contract static when every
    // sequence has an ORF.
    container 'quay.io/biocontainers/cpat:3.0.5--py310h20b60a1_4'

    input:
    tuple val(meta), path(fasta)
    path model
    path hexamer

    output:
    tuple val(meta), path("cpat.ORF_prob.best.tsv"), emit: prob_best
    tuple val(meta), path("cpat.no_ORF.txt")       , emit: no_orf

    script:
    def args = task.ext.args ?: ''
    """
    cpat \\
        -g ${fasta} \\
        -d ${model} \\
        -x ${hexamer} \\
        -o cpat \\
        ${args}

    test -f cpat.no_ORF.txt || touch cpat.no_ORF.txt
    """

    stub:
    """
    touch cpat.ORF_prob.best.tsv cpat.no_ORF.txt
    """
}
