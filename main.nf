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

include { COBS_SEARCH; POSTPROCESS_COBS                               } from './modules/cobs.nf'
include { SKETCH_ASSEMBLY; SKETCH_ANI_DIST; SKETCH_SUBSET; SKETCH_ALL_DIST; SKETCH_TREE; GENERATE_DIST_MATRIX  } from './modules/sketchlib.nf'
include { EXTRACT_ASSEMBLYS_FROM_TAR } from './modules/extract_assembly.nf'
include { PLOT_ANI; PLOT_TREE; SUBSELECT_GRAPH                        } from './modules/plotting.nf'
include { TRIM_TREE                                                   } from './modules/treemmer.nf'

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

    //set up channels for species cobs channel and assembly channels
    channel.fromPath( "${params.cobs_base}/${params.species}*.xz" )
    | set {species_cobs_ch}

    channel.fromPath("${params.assembly_base}/${params.species}*.xz")
    | set { species_assembly_ch }


    manifest = file(params.manifest)

    MANIFEST_PARSE(manifest)
    | combine(species_cobs_ch)
    | COBS_SEARCH
    | groupTuple
    | POSTPROCESS_COBS
    | set { cobs_matches }

    cobs_matches.combine(species_assembly_ch)
    | EXTRACT_ASSEMBLYS_FROM_TAR
    | transpose
    | groupTuple
    | SKETCH_SUBSET
    | set { subset_sketch }
    
    
    SKETCH_ASSEMBLY(MANIFEST_PARSE.out.assemblies)  
    | set { query_sketch }

    subset_sketch.join(query_sketch)
    | SKETCH_ALL_DIST
    | set { all_dists }
       
    subset_sketch.join(query_sketch)
    | SKETCH_ANI_DIST
    | PLOT_ANI

    if (params.sketch_total_ani) {  
        if (params.generate_tree) {             
            SKETCH_TREE(all_dists)
            | PLOT_TREE

            if (params.trim_tree) { 
                TRIM_TREE(SKETCH_TREE.out.tree)
            }
        }

        if (params.cluster_subselection) {
            GENERATE_DIST_MATRIX(all_dists)
            | SUBSELECT_GRAPH
        } 
    }
}
