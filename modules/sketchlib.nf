process SKETCH_ASSEMBLY {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.0.1'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${sketch_db}.skm"), path("${sketch_db}.skd"), emit: assembly_sketch

    script:
    sketch_db = "${meta.ID}_sketch"
    """
    sketchlib sketch -v -k 17 -s 1000 -o ${sketch_db} --threads ${task.cpus} --seq-files ${assembly}
    """
}

process SKETCH_DIST {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.0.1'

    input:
    tuple val(meta), path(query_skm), path(query_skd), path(matches)

    output:
    tuple val(meta), path("dists.txt")

    script:
    """
    sketchlib dist -v -k 17 --ani --threads ${task.cpus} --subset ${matches} ${params.sketchlib_db} > dists.txt
    """

    //sketchlib dist -v -k 17 --ani ${params.sketchlib_db} ${meta.ID}_sketch
}