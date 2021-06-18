////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate input parameters
Workflow.validateWorkflowParams(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )

// Modules: nf-core/modules
include { FASTQC                } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )
include { MULTIQC               } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options              )

// Subworkflows: local
include { INPUT_CHECK           } from '../subworkflows/local/input_check'        addParams( options: [:]                          )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow HICSCAFF {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK (ch_input)
        // - reads : [ meta, reads ]
        // - assemblies : [meta, assembly]

    // Hi-C QC
    FASTQC (INPUT_CHECK.out.reads)
    // Assembly QC
    QUAST(INPUT_CHECK.out.assemblies)           // TODO:: custom
    BUSCO(INPUT_CHECK.out.assemblies)           // TODO:: local
    // Hi-C to Assembly Alignment
    // Alternative - HiC Pro?
    BWAMEM2_INDEX(INPUT_CHECK.out.assemblies)   // TODO:: custom
    BWAMEM2_MEM(INPUT_CHECK.out.reads.join(BWAMEM2_INDEX.out.index))
    // Hi-C QC
    PRESEQ_LCEXTRAP(BWAMEM2_MEM.out.bam)        // What about duplication rate?
    // Hi-C Alignment Contig Orientation processing packages
    // 1) PAIRTOOLS + COOLER
    PAIRTOOLS_PARSE(BWAMEM2_MEM.out.bam)        // Parse and sort
    PAIRTOOLS_DEDUP(PAIRTOOLS_PARSE.out.pairs)  // Deduplicate
    PAIRTOOLS_SELECT(PAIRTOOLS_DEDUP.out.pairs)
    PAIRTOOLS_SPLIT(PAIRTOOLS_SELECT.out.pairs)
    COOLER(PAIRTOOLS_SPLIT.out.pairs)

    // 2) 3D-DNA + JUICER
    BWAMEM2_INDEX()
    JUICER()
    THREEDDNA()
    JBAT()

    // 3) (2) + HiC-Hiker (https://github.com/ryought/hic_hiker)

    // 4) Salsa2
    BEDTOOLS_BAM2BED(BWAMEM2_MEM.out.bam) // Can I sort by name with bam2bed? Maybe samtools sort?
    SAMTOOLS_FAIDX(INPUT_CHECK.out.assemblies)
    SALSA2(BEDTOOLS_BAM2BED.out.bed.join(SAMTOOLS_FAIDX.out.index))


    /*
    Inputs:
        Assemblies
        Hi Libraries - assembly specific? i.e. pairing or crossing input
        Enzyme - library specific? DpnII (dovetail), “Arima”, no enzyme / DNAse (map to juicer options)
        Parameter: 3D-DNA minimum input scaff size
    Processes:
        Index:
        Mapping: special params
        metrics?: assembly in scaffolds > 5Kb, 10Kb, 15Kb
        sitepos script with enzyme (for 3D-DNA input) https://github.com/theaidenlab/juicer/blob/master/misc/generate_site_positions.py
        Dovetail_tools (https://github.com/dovetail-genomics/dovetail_tools)
            preseq lc_extrap which is prone to fail
        Convert pairsam file to formats compatible with 3D-DNA (merged_nodups) and salsa2
        Scaffolding:
            - 3D-DNA
            - salsa2
        Quast
        Busco
    Output:
        .genome file
        Sitepos file
        Stats.txt
        Pairsam.pairs.gz
        .preseq
        Merged_nodups.txt
        Contacts.bed

    */
    BWA_INDEX(
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))


    /*
     * MODULE: Pipeline reporting
     */
    // Get unique list of files containing version information
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS (
        ch_software_versions
    )

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
