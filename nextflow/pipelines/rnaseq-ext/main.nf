#!/usr/bin/env nextflow
/*
 * rnaseq-ext - local extension stages for the adopted nf-core/rnaseq route.
 *
 * Consumes the published results of an adopted nf-core/rnaseq run (see
 * README.md for the contract) and adds the ratified local capabilities:
 * a featureCounts gene-count matrix (with the per-project quantification
 * route switch, decision record D1) and clusterProfiler GO/KEGG/GSEA
 * enrichment on differentialabundance results (D3). Stages are default-off;
 * enable them per project via -params-file.
 */

nextflow.enable.dsl = 2

include { EXT } from './workflows/ext'

workflow {
    EXT()
}
