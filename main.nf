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

include { COBS_SEARCH; POSTPROCESS_COBS                                                              } from './modules/cobs.nf'
include { SKETCH_ASSEMBLY; SKETCH_ANI_DIST; GENERATE_TOTAL_DIST_MATRIX; SKETCH_SUBSET_TOTAL_ANI_DIST } from './modules/sketchlib.nf'
include { BIN_ANI_DISTANCES                                                                          } from './modules/binning.nf'
include { EXTRACT_ASSEMBLYS_FROM_TAR                                                                 } from './modules/extract_assembly.nf'
include { PLOT_ANI; SUBSELECT_GRAPH                                                                  } from './modules/plotting.nf'

//
// SUBWORKFLOWS
//

include { MANIFEST_PARSE } from './subworkflows/manifest_parse.nf'
include { BUILD_TREE     } from './subworkflows/build_tree.nf'

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

    //set up channels for species cobs channel and assembly channels
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
    | set { query_sketch }

    cobs_matches.join(query_sketch)
    | SKETCH_ANI_DIST
    | PLOT_ANI

    BIN_ANI_DISTANCES(SKETCH_ANI_DIST.out.query_ani)
    
    /*
    optional extras
    */

    //using clustering to subselect
    if (params.cluster_subselection) {

        //for this method we need all vs all ANI
        SKETCH_SUBSET_TOTAL_ANI_DIST(cobs_matches)
        | set { subset_ani }

        SKETCH_ANI_DIST.out.query_ani.join(subset_ani)
        | GENERATE_TOTAL_DIST_MATRIX
        | SUBSELECT_GRAPH
    } 
    
    //build a core genome tree for all samples (requires extraction)
    if (params.generate_tree) {
        BUILD_TREE(cobs_matches, query_sketch)
    }
}

