# Blacklist artifact filtering (config key `blacklist`; included only when a
# non-empty path is configured). Filtered copies of the final peak sets (the
# same files annotation consumes) land in 4.peak/blacklist_filtered/ with
# their basenames unchanged; a project-wide before/after count table goes to
# 5.QC/blacklist/. BAM-level filtering is deliberately NOT done (it would
# double the BAM artifacts); peak-level filtering changes FRiP/annotation
# inputs via frip_peak_file()/annot_peak_files() in common.smk.

BLACKLIST_SOURCES = blacklist_sources()


rule blacklist_filter:
    input:
        peaks=lambda wc: BLACKLIST_SOURCES[str(wc)],
        blacklist=BLACKLIST,
    output:
        R("4.peak/blacklist_filtered/{peakfile}"),
    wildcard_constraints:
        peakfile=_group_regex(list(BLACKLIST_SOURCES)),
    log:
        R("logs/blacklist/{peakfile}.log"),
    threads: rthreads("blacklist_filter")
    resources:
        mem_mb=rmem("blacklist_filter"),
        runtime_min=rruntime("blacklist_filter"),
        runtime_sec=rruntime_sec("blacklist_filter"),
    shell:
        """
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} \
            > {output} 2> {log}
        """


rule blacklist_summary:
    input:
        list(BLACKLIST_SOURCES.values())
        + [R(f"4.peak/blacklist_filtered/{b}") for b in BLACKLIST_SOURCES],
    output:
        R("5.QC/blacklist/blacklist_summary.tsv"),
    params:
        srcs=" ".join(list(BLACKLIST_SOURCES.values())),
        dsts=" ".join(R(f"4.peak/blacklist_filtered/{b}") for b in BLACKLIST_SOURCES),
    log:
        R("logs/blacklist/summary.log"),
    threads: rthreads("blacklist_summary")
    resources:
        mem_mb=rmem("blacklist_summary"),
        runtime_min=rruntime("blacklist_summary"),
        runtime_sec=rruntime_sec("blacklist_summary"),
    shell:
        """
        srcs=({params.srcs})
        dsts=({params.dsts})
        printf 'peak_file\\tbefore\\tafter\\tremoved\\n' > {output}
        for i in "${{!srcs[@]}}"; do
            before=$(wc -l < "${{srcs[$i]}}")
            after=$(wc -l < "${{dsts[$i]}}")
            printf '%s\\t%d\\t%d\\t%d\\n' "$(basename "${{dsts[$i]}}")" \\
                "$before" "$after" "$((before - after))" >> {output}
        done
        """
