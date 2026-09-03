# FRiP（Fraction of Reads in Peaks）计算：每 (group, sample) 一份
# ENCODE 参考阈值：TF ≥ 1%（理想 5%+），组蛋白修饰酌情放宽（详见 README）
# 说明：samtools view -c 计数的是 reads（含双端两条），分子分母同口径，比值不变。

FRIP_TSVS = [f"5.QC/frip/{g}__{s}.frip.tsv"
             for g in GROUPS
             for s in GROUPS[g]["treat"] + GROUPS[g]["control"]]

rule frip:
    input:
        bam=lambda wc: sample_bam(wc.sample),
        peaks=lambda wc: group_peak_file(wc.group),
    output:
        "5.QC/frip/{group}__{sample}.frip.tsv",
    wildcard_constraints:
        group=_group_regex(list(GROUPS)),
        sample=_group_regex(SAMPLES),
    log:
        "logs/frip/{group}__{sample}.log",
    threads: 1
    shell:
        """
        mkdir -p 5.QC/frip logs/frip
        total=$(samtools view -c {input.bam})
        inpeak=$(samtools view -c -L {input.peaks} {input.bam})
        awk -v s="{wildcards.sample}" -v g="{wildcards.group}" -v t=$total -v p=$inpeak \
            'BEGIN {{printf "sample\\tgroup\\ttotal_reads\\treads_in_peaks\\tFRiP\\n%s\\t%s\\t%d\\t%d\\t%.4f\\n", s, g, t, p, (t>0 ? p/t : 0)}}' \
            > {output} 2> {log}
        """


rule frip_summary:
    input:
        FRIP_TSVS,
    output:
        tsv="5.QC/frip/FRiP_summary.tsv",
        mqc="5.QC/frip/FRiP_mqc.tsv",
    log:
        "logs/frip/summary.log",
    shell:
        """
        header_written=0
        : > {output.tsv}
        for f in {input}; do
            if [ $header_written -eq 0 ]; then head -n1 "$f" >> {output.tsv}; header_written=1; fi
            tail -n +2 "$f" >> {output.tsv}
        done
        # MultiQC 自定义表（custom content，_mqc.tsv 约定格式）
        {{
            echo "# id: 'frip_table'"
            echo "# section_name: 'FRiP (fraction of reads in peaks)'"
            echo "# format: 'tsv'"
            echo "# plot_type: 'table'"
            echo "# pconfig: {{'id': 'frip_table', 'title': 'FRiP'}}"
            cat {output.tsv}
        }} > {output.mqc} 2> {log}
        """
