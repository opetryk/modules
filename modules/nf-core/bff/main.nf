process BFF {
    publishDir "$params.outdir/$meta.id/$params.mode/hash_demulti/bff", mode:'copy'
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dc475a8f4a175cf1f61e07a020479fbd9b903d1c:e59a2327baad143aef2b1c9283cc605392a7a3ae-0':
        'biocontainers/mulled-v2-dc475a8f4a175cf1f61e07a020479fbd9b903d1c:e59a2327baad143aef2b1c9283cc605392a7a3ae-0' }"

    input:
    tuple val(meta), path(hto_matrix, stageAs: 'hto_data')
    val methods
    val methodsForConsensus
    val cellbarcodeWhitelist
    val metricsFile
    val doTSNE
    val doHeatmap
    val perCellSaturation
    val majorityConsensusThreshold
    val chemistry
    val callerDisagreementThreshold
    val assignmentOutBff
    val preprocess_bff
    val barcodeWhitelist

    output:
    tuple val(meta), path("bff_${meta.id}")

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    mkdir bff_${meta.id}
    """
    template 'bff.R'
}