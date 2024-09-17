//
// Check input samplesheet and get read channels
//

workflow MANIFEST_PARSE {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, "Cannot find path file ${samplesheet}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_assembly_channels(it) }
        .set { assemblies }

    emit:
    assemblies
}

// Function to get list of [ meta, assembly ]
def create_assembly_channels(LinkedHashMap row) {
    def meta = [:]

    meta.ID = row.ID

    def array = []
    // check short reads
    if ( !(row.assembly == 'NA') ) {
        if ( !file(row.assembly).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Assembly file does not exist!\n${row.assembly}"
        }
        assembly = file(row.assembly)
    }

    array = [ meta, assembly ]
    return array
}