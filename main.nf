#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnatoy
========================================================================================
 nf-core/rnatoy Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnatoy
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/sarek -profile docker --input sample.tsv --genome GRCh38"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)

////////////////////////////////////////////////////
/* --        REFERENCES PARAMETER VALUES       -- */
////////////////////////////////////////////////////
/* -- Initialize each params in params.genomes -- */
/* --  catch the command line first if defined -- */
////////////////////////////////////////////////////

params.fasta = Checks.get_genome_attribute(params, 'fasta')
params.gff   = Checks.get_genome_attribute(params, 'gff')

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --             RUN THE WORKFLOW             -- */
////////////////////////////////////////////////////

workflow {

    include { RNATOY } from './workflows/rnatoy' addParams( summary_params: summary_params )
    RNATOY ()

}