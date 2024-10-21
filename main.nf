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
include { LEXICMAP_SEARCH                                                                            } from './modules/lexicmap.nf'
include { SKETCH_ASSEMBLY; SKETCH_ANI_DIST; GENERATE_TOTAL_DIST_MATRIX; SKETCH_SUBSET_TOTAL_ANI_DIST } from './modules/sketchlib.nf'
include { BIN_ANI_DISTANCES                                                                          } from './modules/binning.nf'
include { EXTRACT_ASSEMBLYS_FROM_TAR                                                                 } from './modules/extract_assembly.nf'
include { COLLECT_FILE } from './modules/collect_file.nf'
include { PLOT_ANI; SUBSELECT_GRAPH                                                                  } from './modules/plotting.nf'
include { DOWNLOAD_FASTQS; PUBLISH_FASTQS } from './modules/fastqs.nf'

//from assorted-sub-workflows
include { DOWNLOAD_METADATA                                                                          } from './assorted-sub-workflows/combined_input/modules/ena_downloader.nf'
include { FILTER_METADATA                                                                            } from './assorted-sub-workflows/combined_input/modules/filter_metadata.nf'


//
// SUBWORKFLOWS
//

include { MANIFEST_PARSE } from './subworkflows/manifest_parse.nf'
include { BUILD_TREE     } from './subworkflows/build_tree.nf'
include { QC             } from './assorted-sub-workflows/qc/qc.nf'

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

    /*
    set up channels for species cobs channel and assembly channels
    */

    channel.fromPath( "${params.cobs_base}/${params.species}*.xz" )
    | set {species_cobs_ch}

    channel.fromFilePairs( "${params.sketchlib_db}.{skd,skm}" )
    | collect
    | set { sketchlib_db_ch }

    manifest = file(params.manifest)
    filter_manifest = file(params.filter_manifest, checkIfExists: true)

    //main logic

    MANIFEST_PARSE(manifest)
    | combine(species_cobs_ch)
    | COBS_SEARCH
    | groupTuple
    | POSTPROCESS_COBS
    | set { cobs_matches }


    DOWNLOAD_METADATA(cobs_matches)
    | splitCsv(header: true, sep: "\t")
    | map { meta, full_metadata ->
        def sample_acc = full_metadata.sample_accession
        def cleaned_map = full_metadata.findAll { k, v -> v != '' }
        [ sample_acc, cleaned_map ] //staging sample_acc infront for groupTuple to output from ENADownloader
    }
    | set { sample_metadata }

    FILTER_METADATA(
        DOWNLOAD_METADATA.out.metadata_tsv,
        filter_manifest,
        ["sample_accession"],  // select columns
        true,  // remove_header
        ["sample_accession"] // columns to drop duplicates from
    )
    | set { filtered_cobs_matches }

    SKETCH_ASSEMBLY(MANIFEST_PARSE.out.assemblies)  
    | set { query_sketch }

    SKETCH_ANI_DIST(filtered_cobs_matches.join(query_sketch), sketchlib_db_ch)
     | PLOT_ANI

    BIN_ANI_DISTANCES(SKETCH_ANI_DIST.out.query_ani)
    | splitCsv(header: true, sep: "\t")
    | map { meta, bin_info ->
        def meta_new = [:]
        meta_new.reference_ID = meta.ID //so we can use ID from the sample later
        meta_new.ref_ani_bin = bin_info.ref_ani_bin
        sample = bin_info.query
        
        [ sample, meta_new ] //staging sample infront for groupTuple to output from ENADownloader
    }
    | set{ bin2channel }

    //using clustering to subselect
    if (params.cluster_subselection) {
        //for this method we need all vs all ANI
        bin2channel
        | groupTuple(by: 1)
        | COLLECT_FILE
        | set { samples }

        SKETCH_SUBSET_TOTAL_ANI_DIST(samples, sketchlib_db_ch)
        | GENERATE_TOTAL_DIST_MATRIX
        | SUBSELECT_GRAPH
    }

    bin2channel
    | join(sample_metadata) //replace this with filtered metadata
    | map { join_accession, bin_info, subsampled_metadata ->
        def merged_meta = [:]
        merged_meta = bin_info + subsampled_metadata
        merged_meta.ID = join_accession
        merged_meta
    }
    | filter { it.fastq_ftp.contains(';') } //if its paired its seperated by a semi-colon
    | map{ merged_meta ->
        def (read1_ftp, read2_ftp) = merged_meta.fastq_ftp.split(';')
        def read1_ftp_url = "ftp://${read1_ftp}"
        def read2_ftp_url = "ftp://${read2_ftp}"
        [ merged_meta, read1_ftp_url, read2_ftp_url ]
    }
    | DOWNLOAD_FASTQS
    | set { read_ch }

    read_ch
    | QC
    | filter { it[1] == 'pass' && it[2] == 'pass' }
    | map { it -> it.take(1) } //only keep meta
    | set { filtered_samples }

    filtered_samples
    | join(read_ch)
    | PUBLISH_FASTQS
    
    /*
    optional extras
    */

    if (params.lexicmap_search) {
        channel.fromPath( params.lexicmap_db )
        | set { lexicmap_db_ch }

        LEXICMAP_SEARCH(MANIFEST_PARSE.out.assemblies, lexicmap_db_ch)
    }

    //build a core genome tree for all samples (requires extraction of assemblies)
    if (params.generate_tree) {
        BUILD_TREE(filtered_cobs_matches, query_sketch)
    }
}

