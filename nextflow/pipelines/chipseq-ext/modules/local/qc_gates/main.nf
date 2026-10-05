process QC_GATES_SUMMARY {
    tag "qc_gates"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val thresholds
    val flagstat_suffix
    val dup_suffix
    path flagstats, stageAs: 'flagstat/*'
    path dup_metrics, stageAs: 'dup/*'
    path organelle_summary, stageAs: 'organelle/*'

    output:
    path "qc_gates.tsv"    , emit: tsv
    path "qc_gates_mqc.tsv", emit: mqc

    script:
    def org_max = thresholds.organelle_fraction_max ?: ''
    def dup_arg     = dup_metrics ? "--dup-metrics dup/*" : ''
    def org_arg     = organelle_summary ? "--organelle organelle/*" : ''
    """
    qc_gates_summary.py \\
        --mapping-rate-min ${thresholds.mapping_rate_min} \\
        --dup-rate-max ${thresholds.dup_rate_max} \\
        --organelle-max '${org_max}' \\
        --dup-suffix '${dup_suffix}' \\
        --flagstat-suffix '${flagstat_suffix}' \\
        --out qc_gates.tsv \\
        --mqc qc_gates_mqc.tsv \\
        ${dup_arg} ${org_arg} \\
        flagstat/*
    """

    stub:
    """
    touch qc_gates.tsv qc_gates_mqc.tsv
    """
}
