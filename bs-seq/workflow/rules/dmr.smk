# rules/dmr.smk -- differential methylation (methylKit) over the merged CpG
# tables (v0.2, default off).
#
# Included by workflow/Snakefile after rules/methylation.smk; the rule block
# below is parse-time guarded on config dmr.enabled, so with the default
# configuration the rule is never defined and the DAG is unchanged. One job
# runs workflow/scripts/run_dmr.R over EVERY sample's merged CpG table
# (5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz,
# the bismark 6-column coverage layout that methylKit reads natively via
# pipeline="bismarkCoverage"): pairwise contrasts treat-vs-control against
# dmr.control_group (one output set per contrast, named {treat}_vs_{control}_*),
# plus the tiled DMR table and 6.DMR/DMR_summary.tsv. Config types/ranges and
# the group/batch design are validated at parse time in rules/common.smk
# (validate_config + validate_dmr_design); this rule only exists once those
# checks have passed.
#
# Flag-file pattern (as rna-seq deg): the R script writes every result table,
# the declared output is the 6.DMR/flag.log marker, and the tables themselves
# stay untracked side effects (rerun = delete 6.DMR/ or the flag).

if DMR_ENABLED:

    rule dmr_methylkit:
        input:
            cov=expand(
                R("5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz"),
                sample=SAMPLES,
            ),
            sampleinfo=_SAMPLE_TABLE,
        output:
            flag=R("6.DMR/flag.log"),
        log:
            R("logs/dmr_methylkit.log"),
        params:
            script=os.path.join(WORKFLOW_DIR, "scripts", "run_dmr.R"),
            # Comma-joined (not space-joined) so the R getopt parser receives
            # one argument regardless of the number of files.
            files=lambda wc, input: ",".join(str(p) for p in input.cov),
            control=config["dmr"]["control_group"],
            qvalue=config["dmr"]["qvalue"],
            min_diff=config["dmr"]["min_diff"],
            tile_len=config["dmr"]["tile_len"],
            tile_step=config["dmr"]["tile_step"],
            min_cpg=config["dmr"]["min_cpg"],
            min_cov=config["dmr"]["min_cov"],
            max_cov=config["dmr"]["max_cov"],
            batch=lambda wc: config["dmr"]["batch_correction"],
            outdir=lambda wc, output: os.path.dirname(str(output.flag)),
            rscript=RSCRIPT,
        threads:
            rthreads("dmr")
        resources:
            mem_mb=rmem("dmr"), runtime_min=rruntime("dmr"),
            runtime_sec=rruntime_sec("dmr"),
        shell:
            """
            {params.rscript} {params.script} \\
                -f {params.files} -s {input.sampleinfo} \\
                -r {params.control} -q {params.qvalue} -d {params.min_diff} \\
                -l {params.tile_len} -k {params.tile_step} -g {params.min_cpg} \\
                -c {params.min_cov} -m {params.max_cov} -b {params.batch} \\
                -o {params.outdir} -t {threads} >> {log} 2>&1
            echo `date` " : All done!" > {output.flag}
            """
