/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSentieon.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta not specified!'      }
if (params.known_dbsnp)  { ch_known_dbsnp  = file(params.known_dbsnp)  } else { ch_known_dbsnp  = [] }
if (params.known_mills)  { ch_known_mills  = file(params.known_mills)  } else { ch_known_mills  = [] }
if (params.known_indels) { ch_known_indels = file(params.known_indels) } else { ch_known_indels = [] }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Loaded from modules/local/
//
include { SENTIEON_BWAINDEX } from '../modules/local/sentieon_bwaindex'
include { SENTIEON_BWAMEM   } from '../modules/local/sentieon_bwamem'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_METRICS            } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_LOCUSCOLLECTOR     } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_DEDUP              } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_QUALCAL_RECAL_PRE  } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_QUALCAL_RECAL_POST } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_QUALCAL_RECAL_PLOT } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_HAPLOTYPER         } from '../modules/local/sentieon_driver'
include { SENTIEON_DRIVER as SENTIEON_DRIVER_READWRITER         } from '../modules/local/sentieon_driver'
include { SENTIEON_PLOT   as SENTIEON_PLOT_GCBIAS               } from '../modules/local/sentieon_plot'
include { SENTIEON_PLOT   as SENTIEON_PLOT_QUALDISTRIBUTION     } from '../modules/local/sentieon_plot'
include { SENTIEON_PLOT   as SENTIEON_PLOT_MEANQUALITYBYCYCLE   } from '../modules/local/sentieon_plot'
include { SENTIEON_PLOT   as SENTIEON_PLOT_INSERTSIZEMETRICS    } from '../modules/local/sentieon_plot'
include { SENTIEON_PLOT   as SENTIEON_PLOT_QUALCAL              } from '../modules/local/sentieon_plot'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SENTIEON {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Filter for paired-end reads only
    INPUT_CHECK
        .out
        .reads
        .filter { meta, fasta -> !meta.single_end }
        .set { ch_reads }

    //
    // MODULE: Generate Fasta index file
    //
    SAMTOOLS_FAIDX (
        [ [:], ch_fasta ]
    )
    ch_fai      = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // MODULE: Make BWA index
    //
    SENTIEON_BWAINDEX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(SENTIEON_BWAINDEX.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Sentieon bwa mem command
    //
    SENTIEON_BWAMEM (
        ch_reads,
        ch_fasta,
        ch_fai,
        SENTIEON_BWAINDEX.out.index
    )

    //
    // MODULE: Run Sentieon driver command to get various QC metrics
    //
    SENTIEON_BWAMEM
        .out
        .bam
        .join(SENTIEON_BWAMEM.out.bai)
        .map { it -> it + [ [], [], [], [] ] }
        .set { ch_bam_bai }

    SENTIEON_DRIVER_METRICS (
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        [],
        [],
        []
    )

    //
    // MODULE: Run Sentieon plot to create PDF plots for QC metrics
    //
    SENTIEON_PLOT_GCBIAS (
        SENTIEON_DRIVER_METRICS.out.metrics_gc
    )

    SENTIEON_PLOT_QUALDISTRIBUTION (
        SENTIEON_DRIVER_METRICS.out.metrics_qd
    )

    SENTIEON_PLOT_MEANQUALITYBYCYCLE (
        SENTIEON_DRIVER_METRICS.out.metrics_mq
    )

    SENTIEON_PLOT_INSERTSIZEMETRICS (
        SENTIEON_DRIVER_METRICS.out.metrics_is
    )

    //
    // MODULE: Run Sentieon driver command for LocusCollector
    //
    SENTIEON_DRIVER_LOCUSCOLLECTOR (
        ch_bam_bai,
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(SENTIEON_DRIVER_LOCUSCOLLECTOR.out.versions.first())

    //
    // MODULE: Run Sentieon driver command for Dedup
    //
    ch_bam_bai
        .map { meta, bam, bai, score, score_idx, recal_pre, recal_post -> [ meta, bam, bai ] }
        .join(SENTIEON_DRIVER_LOCUSCOLLECTOR.out.score)
        .join(SENTIEON_DRIVER_LOCUSCOLLECTOR.out.score_idx)
        .map { it -> it + [ [], [] ] }
        .set { ch_bam_bai_score }

    SENTIEON_DRIVER_DEDUP (
        ch_bam_bai_score,
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(SENTIEON_DRIVER_DEDUP.out.versions.first())

    //
    // MODULE: Run Sentieon driver command for QualCal (pre-recalibration)
    //
    SENTIEON_DRIVER_DEDUP
        .out
        .bam
        .join(SENTIEON_DRIVER_DEDUP.out.bai)
        .map { it -> it + [ [], [], [], [] ] }
        .set { ch_dedup_bam_bai }

    SENTIEON_DRIVER_QUALCAL_RECAL_PRE (
        ch_dedup_bam_bai,
        ch_fasta,
        ch_fai,
        ch_known_dbsnp,
        ch_known_mills,
        ch_known_indels
    )
    ch_versions = ch_versions.mix(SENTIEON_DRIVER_QUALCAL_RECAL_PRE.out.versions.first())

    //
    // MODULE: Run Sentieon driver command for QualCal (post-recalibration)
    //
    ch_dedup_bam_bai
        .map { meta, bam, bai, score, score_idx, recal_pre, recal_post -> [ meta, bam, bai, score, score_idx ] }
        .join(SENTIEON_DRIVER_QUALCAL_RECAL_PRE.out.recal_pre)
        .map { it -> it + [ [] ] }
        .set { ch_dedup_recal_bam_bai }

    SENTIEON_DRIVER_QUALCAL_RECAL_POST (
        ch_dedup_recal_bam_bai,
        ch_fasta,
        ch_fai,
        ch_known_dbsnp,
        ch_known_mills,
        ch_known_indels
    )
    ch_versions = ch_versions.mix(SENTIEON_DRIVER_QUALCAL_RECAL_POST.out.versions.first())

    //
    // MODULE: Run Sentieon driver command for QualCal (plot recalibration)
    //
    ch_dedup_recal_bam_bai
        .map { meta, bam, bai, score, score_idx, recal_pre, recal_post -> [ meta, [], [], [], [], recal_pre ] }
        .join(SENTIEON_DRIVER_QUALCAL_RECAL_POST.out.recal_post)
        .set { ch_recal_pre_post }

    SENTIEON_DRIVER_QUALCAL_RECAL_PLOT (
        ch_recal_pre_post,
        [],
        [],
        [],
        [],
        []
    )

    SENTIEON_PLOT_QUALCAL (
        SENTIEON_DRIVER_QUALCAL_RECAL_PLOT.out.recal_csv
    )

    //
    // MODULE: Run Sentieon driver command for Haplotyper
    //
    SENTIEON_DRIVER_HAPLOTYPER (
        ch_dedup_recal_bam_bai,
        ch_fasta,
        ch_fai,
        ch_known_dbsnp,
        [],
        []
    )

    //
    // MODULE: Run Sentieon driver command for ReadWriter
    //
    SENTIEON_DRIVER_READWRITER (
        ch_dedup_recal_bam_bai,
        ch_fasta,
        ch_fai,
        [],
        [],
        []
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSentieon.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
