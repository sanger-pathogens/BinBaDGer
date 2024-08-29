#!/usr/bin/env nextflow

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo

NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


def printHelp() {
    NextflowTool.help_message("", 
                               [],
    params.monochrome_logs, log)
}

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//

include { COBS_SEARCH; POSTPROCESS_COBS } from './modules/cobs.nf'
include { SKETCH_ASSEMBLY; SKETCH_DIST  } from './modules/sketchlib.nf'

//
// SUBWORKFLOWS
//

include { MANIFEST_PARSE   } from './subworkflows/manifest_parse.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    //can decide if we want that logic from the onset
    channel.fromPath( "${params.cobs_base}/${params.species}*.xz" )
    | set {species_cobs_ch}


    manifest = file(params.manifest)

    MANIFEST_PARSE(manifest)
    | combine(species_cobs_ch)
    | COBS_SEARCH
    | groupTuple
    | POSTPROCESS_COBS
    | set { cobs_matches }
    
    SKETCH_ASSEMBLY(MANIFEST_PARSE.out.assemblies)
    | join( cobs_matches )
    | SKETCH_DIST


}
