process TRIM_TREE {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}/trimmed_tree/${meta.ID}", mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/treemmer:a3a1632'

    input:
    tuple val(meta), path(newick)

    output:
    tuple val(meta), path("*.nwk_trimmed_*"), emit: tree

    script:
    """
    Treemmer_v0.3.py ${newick} -X ${params.number_of_leaves} --no_plot
    """
}
