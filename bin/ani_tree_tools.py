#!/usr/bin/env python3

###Code written by John Lees and adapted for purpose here###

import sys
import os
import numpy as np
import pp_sketchlib
import argparse
from tree_builder import generate_phylogeny

import numpy as np

def read_tsv_to_structures(reference_tsv):
    # Initialize the structures
    ref_list = []
    dist_dict = {}

    with open(reference_tsv, 'r') as ref_file:
        for line in ref_file:
            # Split the line by tabs
            sample, reference, ani = line.strip().split('\t')
            if sample not in ref_list:
                ref_list.append(sample)
            if reference not in ref_list:
                ref_list.append(reference)
            # Convert ANI to distance
            dist = 1 - float(ani)
            dist_dict[(sample, reference)] = dist

    # Initialize a square distance matrix
    dist_mat = np.zeros((len(ref_list), len(ref_list)), dtype=np.float32)

    for i, sample in enumerate(ref_list):
        for j in range(i, len(ref_list)):
            reference = ref_list[j]
            if i == j:
                dist = 0.0  # Distance to self is 0
            else:
                dist_mat[i, j] = dist_dict[(sample, reference)]

    return ref_list, dist_mat

def read_tsv_to_core_accession(reference_tsv):
    # Initialize the structures
    ref_list = []
    distances = []

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
            distances.append((sample, reference, float(core_dist), float(acc_dist)))

    # Initialize matrices
    num_samples = len(ref_list)
    core_dist_mat = np.zeros((num_samples, num_samples))
    acc_dist_mat = np.zeros((num_samples, num_samples))

    # Update the distance matrices
    for sample, reference, core_dist, acc_dist in distances:
        sample_idx = ref_list.index(sample)
        reference_idx = ref_list.index(reference)

        core_dist_mat[sample_idx, reference_idx] = core_dist
        core_dist_mat[reference_idx, sample_idx] = core_dist

        acc_dist_mat[sample_idx, reference_idx] = acc_dist
        acc_dist_mat[reference_idx, sample_idx] = acc_dist

    return ref_list, core_dist_mat, acc_dist_mat

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
    parser.add_argument('--core_accession', action='store_true', help='parse input TSV as core + accession rather than single ANI scores')
    args = parser.parse_args()

    if args.phylip_path:
        phylip_path = args.phylip_path
        sys.stderr.write(f"Using provided PHYLIP file: {phylip_path}\n")
    else:
        # Read the TSV and process data
        if args.core_accession:
            ref_list, dist_mat, acc_dist_mat = read_tsv_to_core_accession(args.dist_tsv_path)        
        else:
            ref_list, dist_mat = read_tsv_to_structures(args.dist_tsv_path)

        # Generate the PHYLIP matrix
        phylip_path = generate_phylip_matrix(ref_list, dist_mat, args.meta_ID)


    # Conditionally build the tree if --build_tree is specified
    if args.build_tree:
        generate_phylogeny(phylip_path, args.meta_ID, "nwk", True)
    else:
        sys.stderr.write("Skipping tree generation\n")
