#!/usr/bin/env nextflow
/*
 * chipseq-ext - local extension stages for the adopted nf-core/chipseq route.
 *
 * Consumes the published results of an adopted nf-core/chipseq run (see
 * README.md for the contract) and adds the ratified local capabilities:
 * organelle QC, QC gate table, spike-in scaling, narrow-peak IDR, pooled
 * SEACR calling, HOMER motif discovery, and DiffBind differential binding.
 * All stages are default-off; enable them per project via -params-file.
 */

nextflow.enable.dsl = 2

include { EXT } from './workflows/ext'

workflow {
    EXT()
}
