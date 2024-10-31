process COBS_SEARCH {
    tag "${meta.ID}"
    label "cpu_2"
    label "mem_10"
    label "time_30m"

    container 'quay.io/biocontainers/cobs:0.3.0--hdcf5f25_1'

    input:
    tuple val(meta), path(assembly), path(index)

    output:
    tuple val(meta), path("${meta.ID}_${index}_matches.txt"), emit: matched_genomes

    script:
    //I want to stream - Calculate uncompressed size with xz before with --list and use that to allocate memory for cobs to stream from the index
    """
    uncompressed_size=\$(xz --robot --list -vv ${index} | awk '\$1 == "totals" {print \$5}')

    cobs query --load-complete -T ${task.cpus} --threshold ${params.cobs_threshold} -i <(xzcat --no-sparse --ignore-check ${index}) --index-sizes \$uncompressed_size -f ${assembly} > ${meta.ID}_${index}_matches.txt
    """
}

process POSTPROCESS_COBS {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), path(matches)

    output:
    tuple val(meta), path("combined_matches.txt"), emit: combined_cobs_matches

    script:
    """
    postprocess_cobs.py --selection_method ${params.selection_method} --seed 1234 -n ${params.number_of_cobs_matches} ${matches}
    """
}
