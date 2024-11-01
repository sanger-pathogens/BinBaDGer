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

    publishDir "${params.outdir}/tree/${meta.ID}", pattern: '*.png', mode: 'copy', overwrite: true

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
    tag "${meta.reference_ID}_${meta.ref_ani_bin}"
    label "cpu_1"
    label "mem_8"
    label "time_12"

    publishDir "${params.outdir}/clusters/${meta.reference_ID}/${meta.ref_ani_bin}", pattern: "*.png", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/clusters/${meta.reference_ID}/${meta.ref_ani_bin}", pattern: "*.gif", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/clusters/${meta.reference_ID}/${meta.ref_ani_bin}", pattern: "*.txt", mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/python_graphics:1.1.4'

    input:
    tuple val(meta), path(phylip)

    output:
    path("*.txt"), emit: representatives, optional: true
    path("*.csv"), emit: clusters, optional: true
    path("*.png"), optional: true
    path("*.gif"), optional: true

    script:
    def make_gif = params.make_gif ? "--plot_selection_plots" : ""
    def representatives = params.representatives ? "--n_representatives ${params.representatives}" : ""
    """
    subselect_graph.py --phylip ${phylip} --methods ${params.cluster_method} ${make_gif} ${representatives}
    """
}
