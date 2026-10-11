process FINAL_LNCRNA_SET {
    tag "final-set"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(pfam_ids, stageAs: 'pfam/*')
    tuple val(meta2), path(diamond_ids, stageAs: 'diamond/*')
    tuple val(meta3), path(noncoding_gtf)
    path reference_gtf
    path genome

    output:
    tuple val(meta), path("final_lncRNA.fa")          , emit: final_fasta
    tuple val(meta), path("final_lncRNA.gtf")         , emit: final_gtf
    tuple val(meta), path("ref_withLncRNA.gtf")       , emit: ref_with_lncrna_gtf
    tuple val(meta), path("homology.hitIDs.txt")      , emit: hit_ids

    script:
    def pfam_arg   = pfam_ids ? "--pfam 'pfam/*'" : ''
    def diamond_arg = diamond_ids ? "--diamond 'diamond/*'" : ''
    """
    final_lncrna.py \\
        --gtf ${noncoding_gtf} \\
        --reference ${reference_gtf} \\
        --genome ${genome} \\
        ${pfam_arg} ${diamond_arg} \\
        --out-gtf final_lncRNA.gtf \\
        --out-fa final_lncRNA.fa \\
        --out-ref ref_withLncRNA.gtf \\
        --out-hits homology.hitIDs.txt
    """

    stub:
    """
    touch final_lncRNA.fa final_lncRNA.gtf ref_withLncRNA.gtf homology.hitIDs.txt
    """
}
