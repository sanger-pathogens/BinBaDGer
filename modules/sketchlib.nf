process SKETCH_ASSEMBLY {
    tag "${meta.ID}"
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.2_sd28_fix'

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
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.2_sd28_fix'

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
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.2_sd28_fix'

    input:
    // the query_skm and the query_skd relate to parts of a sketchlib database. The skm is the metadata
    // and the skd is the data itself. These together make a sketch which can be compared to other sketches
    // or, in the case of a multisketch, pairwise across all assemblies or reads within that sketch.
    // https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/indexes/README.md
    tuple val(meta), path(subset), path(query_skm), path(query_skd)
    tuple val(sketchlib_db), path(sketchlib_db_files) //bring them along

    output:
    tuple val(meta), path("${meta.ID}_ani_data.tsv"), emit: query_ani

    script:
    query_db = "${meta.ID}_sketch"
    """
    sketchlib dist -v -k 17 --subset ${subset} --ani ${sketchlib_db}.skm ${query_db} > ${meta.ID}_ani_data.tsv
    """
}

process SKETCH_SUBSET_TOTAL_ANI_DIST {
    tag "${meta.reference_ID}_${meta.ref_ani_bin}"
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.2_sd28_fix'

    input:
    tuple val(meta), path(subset)
    tuple val(sketchlib_db), path(sketchlib_db_files)

    output:
    tuple val(meta), path("${meta.reference_ID}_${meta.ref_ani_bin}_betweenness_ani.tsv"), emit: subset_ani

    script:
    """
    sketchlib dist -v -k 17 --subset ${subset} --ani ${sketchlib_db}.skm > ${meta.reference_ID}_${meta.ref_ani_bin}_betweenness_ani.tsv
    """
}

process SKETCH_CORE_ACC_DIST {
    tag "${meta.ID}"
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/ssd28/experimental/pp-sketchlib-rust:0.1.2_sd28_fix'

    input:
    tuple val(meta), path(ref_skm), path(ref_skd), path(query_skm), path(query_skd)

    output:
    tuple val(meta), path("${meta.ID}_total_core_data.tsv"), emit: total_ani

    script:
    query_db = "${meta.ID}_sketch"
    """
    sketchlib dist -v ${ref_skm} ${query_db} > ${meta.ID}_total_core_data.tsv
    sketchlib dist -v ${ref_skm} >> ${meta.ID}_total_core_data.tsv
    """
}

process SKETCH_TREE {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_8"
    label "time_1"

    //this is all testing so using experimental repo
    container 'quay.io/ssd28/experimental/rapidnj:2.3.2-c1'

    publishDir "${params.outdir}/tree/${meta.ID}", pattern: '*.nwk', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(total_tsv)

    output:
    tuple val(meta), path("*.nwk"), emit: tree

    script:
    """
    ani_tree_tools.py --dist_tsv_path ${total_tsv} --meta_ID ${meta.ID} --core_accession --build_tree
    """
}

process GENERATE_TOTAL_DIST_MATRIX {
    tag "${meta.reference_ID}_${meta.ref_ani_bin}"
    label "cpu_4"
    label "mem_8"
    label "time_1"

    container 'quay.io/ssd28/experimental/rapidnj:2.3.2-c1'

    input:
    tuple val(meta), path(betweenness_tsv)

    output:
    tuple val(meta), path("*.phylip"), emit: matrix

    script:
    """
    ani_tree_tools.py --dist_tsv_path ${betweenness_tsv} --meta_ID ${meta.ref_ani_bin}
    """
}
