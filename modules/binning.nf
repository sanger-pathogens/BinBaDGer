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
    path("${meta.ID}_binning.log")

    script:
    retain_below_bins = params.retain_below_bins ? '--assign_outsiders' : ''
    """
    bin_ani.py --input_tsv ${ani_tsv} --output_tsv ${meta.ID}_binned.tsv --log_file ${meta.ID}_binning.log --bins ${params.bin_ranges} ${retain_below_bins}
    """
}