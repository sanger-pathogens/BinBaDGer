process EXTRACT_ASSEMBLYS_FROM_TAR {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"


    container 'quay.io/sangerpathogens/python_graphics:1.1.3'

    input:
    tuple val(meta), path(matches), path(tar_file)

    output:
    tuple val(meta), path("*.fasta"), emit: fastas, optional: true

    script:
    """
    tar_extract.py -f ${tar_file} -n ${matches}
    """
}