#!/usr/bin/env nextflow
/*
 * atacseq-ext - local extension stages for the adopted nf-core/atacseq route.
 *
 * Consumes the published results of an adopted nf-core/atacseq run (see
 * README.md for the contract) and adds the ratified local capabilities:
 * organelle QC, QC gate table, spike-in scaling, and the per-group TOBIAS
 * footprinting chain. All stages are default-off; enable them per project
 * via -params-file.
 */

nextflow.enable.dsl = 2

include { EXT } from './workflows/ext'

workflow {
    EXT()
}
