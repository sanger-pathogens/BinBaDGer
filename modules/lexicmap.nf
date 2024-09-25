process LEXICMAP_SEARCH {
    tag "${meta.ID}"
    cpus = 48 //yes really
    label "mem_96"
    label "time_12"

    container 'quay.io/sangerpathogens/lexicmap:0.4.0'

    publishDir "${params.outdir}/lexicmap/${meta.ID}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.ID}_lexicmap_matches.tsv"), emit: matched_genomes

    script:
    """
    lexicmap search -d ${params.lexicmap_db} ${assembly} \\
        -o ${meta.ID}_lexicmap_matches.tsv \\
        --min-qcov-per-hsp 70 \\
        --min-qcov-per-genome 70 \\
        --top-n-genomes 1000
    """
}