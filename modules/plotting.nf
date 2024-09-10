process PLOT_ANI {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}/plots/${meta.ID}", mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/python_graphics:1.1.3'

    input:
    tuple val(meta), path(ani_tsv)

    output:
    tuple val(meta), path("*.png"), emit: plots

    script:
    """
    plot_ani.py ${ani_tsv}
    """
}

process PLOT_TREE {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}/tree/${meta.ID}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/experimental/rapidnj:2.3.2-c1'

    input:
    tuple val(meta), path(newick)

    output:
    tuple val(meta), path("*.png"), emit: plots

    script:
    """
    plot_tree.py ${newick} ${meta.ID}_tree.png
    """
}

process SUBSELECT_GRAPH {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}/clusters/${meta.ID}", mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/python_graphics:1.1.3'

    input:
    tuple val(meta), path(phylip)

    output:
    tuple val(meta), path("*.csv"), emit: clusters
    tuple val(meta), path("*.png")
    tuple val(meta), path("*.txt")

    script:
    """
    subselect_graph.py --phylip ${phylip} --methods all > log.txt
    """
}