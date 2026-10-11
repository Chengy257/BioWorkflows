#!/usr/bin/env nextflow
/*
 * lncrna-nf - BioWorkflows lncRNA identification workflow (local DSL2).
 *
 * Consumes the adopted nf-core/rnaseq route's published BAMs (the Stage-2
 * handoff) and runs the ratified identification chain: StringTie assembly
 * -> merge against the protein-coding GTF -> gffcompare classcode-u novel
 * extraction -> coding potential (CPC2 AND CPAT, length >= 200) -> Pfam +
 * DIAMOND/NR homology filters -> final lncRNA set -> featureCounts against
 * ref_withLncRNA.gtf. Scientific record:
 * nextflow/docs/lncrna-scientific-decision-record.md.
 */

nextflow.enable.dsl = 2

include { IDENTIFY } from './workflows/identify'

workflow {
    IDENTIFY()
}
