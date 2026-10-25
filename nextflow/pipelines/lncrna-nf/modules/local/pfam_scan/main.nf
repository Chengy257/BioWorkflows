process PFAM_SCAN {
    tag "pfam-scan"
    label 'process_medium'
    // pfam_scan.pl is not packaged as an nf-core module; biocontainers
    // carries pfam_scan 1.6 (verified on quay 2026-10-12). The 1.6 image
    // lacks the standalone 'translate' binary that -translate needs
    // (found at the 2026-10-26 real-container round), so the nucleotide
    // FASTA is translated here (six-frame ORFs) and pfam_scan consumes
    // protein sequences; hit ids carry the '<tx>|orfN' separator the hit
    // parser strips.
    container 'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_5'

    input:
    tuple val(meta), path(novel_gtf)
    tuple val(meta2), path(fasta)
    path db_dir, stageAs: 'pfam_data'
    val evalue

    output:
    tuple val(meta), path("pfam_scan.output"), emit: output
    tuple val(meta), path("pfam.hitIDs.txt") , emit: hit_ids

    script:
    """
    six_frame_orfs.py ${fasta} orfs.fa --min-aa 20
    pfam_scan.pl \\
        -fasta orfs.fa \\
        -dir pfam_data \\
        -as \\
        -cpu ${task.cpus} \\
        -e_seq ${evalue} \\
        > pfam_scan.output
    pfam_hits.py \\
        --pfam pfam_scan.output \\
        --evalue ${evalue} \\
        --out pfam.hitIDs.txt
    """

    stub:
    """
    touch pfam_scan.output pfam.hitIDs.txt
    """
}
