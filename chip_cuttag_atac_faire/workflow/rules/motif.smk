# HOMER motif enrichment on the final peak sets (v0.5; included only when
# motif.enabled is true). De novo + known motif discovery per group via
# findMotifsGenome.pl; HOMER is an external distribution (configureHomer)
# resolved through the software.yaml paths: mechanism (CHIP_HOMER_FINDMOTIFS)
# or left as findMotifsGenome.pl on PATH. motif.homer_genome names the HOMER
# genome tag (e.g. hg38, mm10) or a custom:/path/to/genome directory and is
# mandatory when the stage is enabled. Runs can take hours — see the
# motif_enrichment resource entry before queueing large peak sets.

rule motif_enrichment:
    input:
        peaks=lambda wc: motif_peak_file(wc.group),
    output:
        directory(R("6.motif/{group}")),
    wildcard_constraints:
        group=_group_regex(list(GROUPS)),
    params:
        homer=HOMER_BIN,
        genome=MOTIF["homer_genome"],
        size=MOTIF["size"],
        background="" if not MOTIF["background"] else f"-bg {MOTIF['background']}",
        extra=MOTIF["extra"],
    log:
        R("logs/motif/{group}.log"),
    threads: rthreads("motif_enrichment")
    resources:
        mem_mb=rmem("motif_enrichment"),
        runtime_min=rruntime("motif_enrichment"),
        runtime_sec=rruntime_sec("motif_enrichment"),
    shell:
        """
        mkdir -p {output}
        # findMotifsGenome.pl calls its helper tools (homer2, findKnownMotifs.pl,
        # compareMotifs.pl, ...) by bare name, and HOMER continues "successfully"
        # with partial output when they are not on PATH — prepend the HOMER bin
        # directory (dirname of the resolved entry point) so the full pipeline
        # actually runs (WSL real-run finding 2026-09-09)
        homer_bindir="$(dirname "{params.homer}")"
        if [ "$homer_bindir" != "." ]; then export PATH="$homer_bindir:$PATH"; fi
        {params.homer} {input.peaks} {params.genome} {output} \
            -size {params.size} {params.background} -cpu {threads} \
            {params.extra} > {log} 2>&1
        """
