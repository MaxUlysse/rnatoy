// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process TOPHAT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::tophat=2.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tophat:2.1.1--py27_3"
    } else {
        container "quay.io/biocontainers/tophat:2.1.1--py27_3"
    }

    input:
    path genome
    path gff
    path index
    tuple val(meta), path(reads)
 
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"          , emit: version

    script:
    def software  = getSoftwareName(task.process)

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    ln -s $genome bowtie2/.

    tophat2 -p ${task.cpus} --GTF $gff \$INDEX $reads
    mv tophat_out/accepted_hits.bam .

    echo "2.1.1" > ${software}.version.txt
    """
}