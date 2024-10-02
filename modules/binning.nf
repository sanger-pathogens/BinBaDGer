process BIN_ANI_DISTANCES {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container 'quay.io/sangerpathogens/python_graphics:1.1.3'

    publishDir "${params.outdir}/bins/${meta.ID}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(ani_tsv)

    output:
    tuple val(meta), path("${meta.ID}_binned.tsv"), emit: binned_ani

    script:
    """
    bin_ani.py ${ani_tsv} ${meta.ID}_binned.tsv 3
    """
}