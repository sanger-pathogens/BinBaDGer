process COBS_SEARCH {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_16"
    label "time_1"

    container 'quay.io/biocontainers/cobs:0.3.0--hdcf5f25_1'

    input:
    tuple val(meta), path(assembly), path(index)

    output:
    tuple val(meta), path("${meta.ID}_${index}_matches.txt"), emit: matched_genomes

    script:
    //I want to stream - Calculate uncompressed size with xz before with --list and use that to allocate memory for cobs to stream from the index 
    """
    uncompressed_size=\$(xz --robot --list -vv ${index} | awk '\$1 == "totals" {print \$5}')

    cobs query --load-complete -T ${task.cpus} -i <(xzcat --no-sparse --ignore-check ${index}) --index-sizes \$uncompressed_size -f ${assembly} > ${meta.ID}_${index}_matches.txt
    """
}

process POSTPROCESS_COBS {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), path(matches)

    output:

    script:
    //I want to stream - Calculate uncompressed size with xz before with --list and use that to allocate memory for cobs to stream from the index 
    """
    postprocess_cobs.py -n ${params.samples} ${matches} > combined_matches.txt
    """
}