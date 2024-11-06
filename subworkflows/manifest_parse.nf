//
// Check input samplesheet and get read channels
//

workflow MANIFEST_PARSE {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    def seenIds = new HashSet()

    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, "Cannot find path file ${samplesheet}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_assembly_channels(it, seenIds) }
        .set { assemblies }

    emit:
    assemblies
}

// Function to get list of [ meta, assembly ]
def create_assembly_channels(LinkedHashMap row, HashSet seenIds) {
    def meta = [:]

    meta.ID = row.ID

    // Check if ID already exists
    if (!seenIds.add(row.ID)) {
        error("ERROR: Duplicate ID found in samplesheet! '${row.ID}' is used multiple times.")
    }

    def assembly
    if ( !file(row.assembly).exists() ) {
        error("ERROR: Please check input samplesheet -> Assembly file does not exist!\n${row.assembly}")
    }
    assembly = file(row.assembly)

    return [ meta, assembly ]
}