process DOWNLOAD_FASTQS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    maxForks 10

    publishDir "${params.outdir}/${meta.ID}/fastqs", mode: 'copy', overwrite: true, enabled: params.output_all_fastqs

    input:
    tuple val(meta), val(fastq_path_1), val(fastq_path_2)

    output:
    tuple val(meta), path(read_1), path(read_2), emit: fastqs

    script:
    read_1 = "${meta.ID}_1.fastq.gz"
    read_2 = "${meta.ID}_2.fastq.gz"
    """
    wget --progress=dot:giga \\
            -O ${read_1} \\
            ${fastq_path_1}

    wget --progress=dot:giga \\
            -O ${read_2} \\
            ${fastq_path_2}
    """
}

process PUBLISH_FASTQS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    publishDir "${params.outdir}/${meta.ID}/fastqs", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path(read_1), path(read_2), emit: fastqs

    script:
    """
    """
}
