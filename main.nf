#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/sentieon
========================================================================================
    Github : https://github.com/nf-core/sentieon
    Website: https://nf-co.re/sentieon
    Slack  : https://nfcore.slack.com/channels/sentieon
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta        = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fai          = WorkflowMain.getGenomeAttribute(params, 'fai')
params.bwa_index    = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.known_dbsnp  = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp')
params.known_mills  = WorkflowMain.getGenomeAttribute(params, 'known_mills')
params.known_indels = WorkflowMain.getGenomeAttribute(params, 'known_indels')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SENTIEON } from './workflows/sentieon'

//
// WORKFLOW: Run main nf-core/sentieon analysis pipeline
//
workflow NFCORE_SENTIEON {
    SENTIEON ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SENTIEON ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
