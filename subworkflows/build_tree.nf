/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//

include { EXTRACT_ASSEMBLYS_FROM_TAR                       } from '../modules/extract_assembly.nf'
include { SKETCH_SUBSET; SKETCH_CORE_ACC_DIST; SKETCH_TREE } from '../modules/sketchlib.nf'
include { PLOT_TREE                                        } from '../modules/plotting.nf'
include { TRIM_TREE                                        } from '../modules/treemmer.nf'


workflow BUILD_TREE {
    take:
    cobs_matches
    query_sketch

    main:
    
    channel.fromPath("${params.assembly_base}/${params.species}*.xz")
    | set { species_assembly_ch }

    cobs_matches.combine(species_assembly_ch)
    | EXTRACT_ASSEMBLYS_FROM_TAR
    | transpose
    | groupTuple
    | SKETCH_SUBSET
    | set { subset_sketch }

    subset_sketch.join(query_sketch)
    | SKETCH_CORE_ACC_DIST     
    | SKETCH_TREE
    | PLOT_TREE

    if (params.trim_tree) { 
        TRIM_TREE(SKETCH_TREE.out.tree)
    }

}