process PREPROCESS {

    publishDir "${params.outdir}/preprocess", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.csv"), emit: readnum

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    concatenate.R samplesheet nanopore ${prefix}.csv
    """
}

