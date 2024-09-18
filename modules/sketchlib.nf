process SKETCH_ASSEMBLY {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.1_sd28_fix'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${sketch_db}.skm"), path("${sketch_db}.skd"), emit: assembly_sketch

    script:
    sketch_db = "${meta.ID}_sketch"
    """
    sketchlib sketch -v -k 3,17,35 -o ${sketch_db} -s 1024 --seq-files ${assembly}
    """
}

process SKETCH_SUBSET {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.1_sd28_fix'

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("${sketch_db}.skm"), path("${sketch_db}.skd"), emit: assemblies_sketch

    script:
    sketch_db = "${meta.ID}_sketch_subset"
    """
    sketchlib sketch -v -k 3,17,35 -o ${sketch_db} -s 1024 --seq-files ${assemblies}
    """
}

process SKETCH_ANI_DIST {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.1_sd28_fix'

    input:
    tuple val(meta), path(ref_skm), path(ref_skd)
    tuple val(meta), path(query_skm), path(query_skd)

    output:
    tuple val(meta), path("${meta.ID}_ani_data.tsv"), emit: query_ani

    script:
    query_db = "${meta.ID}_sketch"
    """
    sketchlib dist -v -k 17 --ani ${ref_skm} ${query_db} > ${meta.ID}_ani_data.tsv
    """
}

process SKETCH_ALL_DIST {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.1_sd28_fix'

    input:
    tuple val(meta), path(ref_skm), path(ref_skd), path(query_skm), path(query_skd)

    output:
    tuple val(meta), path("${meta.ID}_total_core_data.tsv"), emit: total_ani

    script:
    query_db = "${meta.ID}_sketch"
    """
    sketchlib dist -v ${ref_skm} ${query_db} > ${meta.ID}_query_core_data.tsv
    sketchlib dist -v ${ref_skm} > ${meta.ID}_ref_core_data.tsv
    
    cat ${meta.ID}_ref_core_data.tsv ${meta.ID}_query_core_data.tsv > ${meta.ID}_total_core_data.tsv 
    """
}

process SKETCH_TREE {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_8"
    label "time_1"

    //this is all testing so using experimental repo
    container 'quay.io/ssd28/experimental/rapidnj:2.3.2'

    publishDir "${params.outdir}/tree/${meta.ID}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(total_tsv)

    output:
    tuple val(meta), path("*.nwk"), emit: tree

    script:
    """
    ani_tree_tools.py --dist_tsv_path ${total_tsv} --meta_ID ${meta.ID} --build_tree
    """
}

process GENERATE_DIST_MATRIX {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_8"
    label "time_1"

    //for ppsketchlib
    container 'quay.io/ssd28/experimental/rapidnj:2.3.2-c1'

    publishDir "${params.outdir}/tree/${meta.ID}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(total_tsv)

    output:
    tuple val(meta), path("*.phylip"), emit: matrix

    script:
    """
    ani_tree_tools.py --dist_tsv_path ${total_tsv} --meta_ID ${meta.ID}
    """
}