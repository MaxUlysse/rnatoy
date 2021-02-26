////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

include { extract_fastq } from '../modules/local/process/functions.nf'

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.gff
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
if (params.gff)   { ch_gff   = file(params.gff) }   else { exit 1, 'Genome gff file not specified!' }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

def fastqc_options         = modules['fastqc']

include { FASTQC                             } from '../modules/nf-core/software/fastqc/main.nf'        addParams( options: fastqc_options )
include { BOWTIE2_BUILD                      } from '../modules/nf-core/software/bowtie2/build/main.nf' addParams( options: multiqc_options )
include { MULTIQC                            } from '../modules/nf-core/software/multiqc/main.nf'       addParams( options: multiqc_options )
include { TOPHAT                             } from '../modules/local/process/tophat'                   addParams( options: [:] )
include { GET_SOFTWARE_VERSIONS              } from '../modules/local/process/get_software_versions'    addParams( options: [publish_files : ['csv':'']] )
include { INPUT_CHECK                        } from '../modules/local/subworkflow/input_check'          addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

workflow RNATOY {

    ch_software_versions = Channel.empty()

    INPUT_CHECK ( 
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_fastq }

    FASTQC( ch_fastq )

    BOWTIE2_BUILD ( ch_fasta )

    TOPHAT( ch_fasta, ch_gff, BOWTIE2_BUILD.out.index, ch_fastq )

    ch_software_versions = ch_software_versions.mix(
        FASTQC.out.version.first().ifEmpty(null),
        BOWTIE2_BUILD.out.version.ifEmpty(null),
        TOPHAT.out.version.ifEmpty(null)
    )

    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        GET_SOFTWARE_VERSIONS.out.yaml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        FASTQC.out.zip.collect{it[1]}.ifEmpty([])
    )
}
