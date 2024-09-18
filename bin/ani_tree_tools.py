#!/usr/bin/env python3

###Code written by John Lees and adapted for purpose here###

import sys
import os
import numpy as np
import pp_sketchlib
import argparse
from tree_builder import generate_phylogeny

def read_tsv_to_structures(reference_tsv):
    # Initialize the structures
    ref_list = []
    dist_mat = []

    # Read and process the reference TSV
    with open(reference_tsv, 'r') as ref_file:
        for line in ref_file:
            # Split the line by tabs
            sample, reference, core_dist, acc_dist = line.strip().split('\t')
            # Add the unique samples to ref_list
            if sample not in ref_list:
                ref_list.append(sample)
            if reference not in ref_list:
                ref_list.append(reference)
            # Convert ANI to distances and add to dist_mat

            dist_mat.append([float(core_dist), float(core_dist)])
            

    # Convert dist_mat to a numpy array
    dist_mat = np.array(dist_mat, dtype=np.float32)


    return ref_list, dist_mat

def update_distance_matrices(dist_mat, threads = 1):
    """Convert distances from long form (1 matrix with n_comparisons rows and 2 columns)
    to a square form (2 NxN matrices)
    Args:
        ref_list (list)
            List of references
        dist_mat (numpy.array)
            Two column long form list of core and accessory distances
            for pairwise comparisons between reference db sequences
        threads (int)
            Number of threads to use
    Returns:
        seqLabels (list)
            Combined list of reference and query sequences
        matrix (numpy.array)
            NxN array of core distances for N sequences
    """
    matrix = pp_sketchlib.longToSquare(dist_mat[:, [0]], threads)
    return matrix

def generate_phylip_matrix(ref_list, matrix, meta_ID):
    # generate phylip matrix
    phylip_name = f"{meta_ID}_distances.phylip"
    with open(phylip_name, 'w') as pFile:
        pFile.write(str(len(ref_list)) + "\n")
        for Dist, ref in zip(matrix, ref_list):
            pFile.write(ref)
            pFile.write(' ' + ' '.join(map(str, Dist)))
            pFile.write("\n")
    full_path = os.path.abspath(phylip_name)
    return full_path
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process input files")
    parser.add_argument('-r', '--dist_tsv_path', type=str, required=True, help='Input TSV file with Reference pairwise ANI data')
    parser.add_argument('--meta_ID', type=str, required=True, help='ID of dataset')
    parser.add_argument('--build_tree', action='store_true', help='Option to build tree')
    parser.add_argument('--phylip_path', type=str, help='Optional: Pre-generated PHYLIP file path')
    args = parser.parse_args()

    if args.phylip_path:
        phylip_file = args.phylip_path
        sys.stderr.write(f"Using provided PHYLIP file: {phylip_file}\n")
    else:
        # Read the TSV and process data
        ref_list, dist_mat = read_tsv_to_structures(args.dist_tsv_path)
        matrix = update_distance_matrices(dist_mat)

        # Generate the PHYLIP matrix
        phylip_path = generate_phylip_matrix(ref_list, matrix, args.meta_ID)

    # Conditionally build the tree if --build_tree is specified
    if args.build_tree:
        generate_phylogeny(phylip_path, args.meta_ID, "nwk", True)
    else:
        sys.stderr.write("Skipping tree generation\n")