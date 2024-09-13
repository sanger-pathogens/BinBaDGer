#!/usr/bin/env python3

import tarfile
import logging
from pathlib import Path
import argparse

# Function to iterate over a batch of TAR files and extract relevant FASTA files
def iterate_over_batch(asms_fn, selected_rnames):
    """Iterate over an xz-compressed TAR file corresponding to a batch with individual FASTA files.

    Args:
        asms_fn (str): xz-compressed TAR file with FASTA files.
        selected_rnames (list): Set of selected FASTA files for which a Minimap instance will be created (note: can contain rnames from other batches).
    Returns:
        (rname (str), rfa (str))
    """
    logging.info(f"Opening {asms_fn}")
    skipped = 0

    with tarfile.open(asms_fn, mode="r:xz") as tar:
        for member in tar.getmembers():
            # extract file headers
            name = member.name
            rname = Path(name).stem
            if rname not in selected_rnames:
                logging.debug(f"Skipping {rname} ({name})")
                skipped += 1
                continue
            # extract file content
            if skipped > 0:
                logging.info(f"Skipping {skipped} references in {asms_fn} (no hits for them)")
                skipped = 0
            logging.info(f"Extracting {rname} ({name})")
            f = tar.extractfile(member)
            rfa = f.read()
            yield rname, rfa

    if skipped > 0:
        logging.info(f"Skipping {skipped} references in {asms_fn}")

# Function to save FASTA content to a file
def save_fasta(rname, rfa):
    """Save the extracted sequence to a file named after the accession (rname)."""
    filename = f"{rname}.fasta"
    logging.info(f"Saving {rname} to {filename}")
    with open(filename, "wb") as fasta_file:
        fasta_file.write(rfa)

# Helper function to extract sequences from multiple files based on index and save them
def extract_sequences_from_files(file_list, selected_rnames):
    for tar in file_list:
        for rname, rfa in iterate_over_batch(tar, selected_rnames):
            save_fasta(rname, rfa)

# Function to read lists from files
def read_list_from_file(filepath):
    """Read lines from a file, stripping newlines."""
    with open(filepath, "r") as f:
        return [line.strip() for line in f.readlines()]

# Argument parser for the script
def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract FASTA sequences from xz-compressed TAR files.")

    parser.add_argument(
        "-f", "--files",
        nargs="+",  # Accept multiple tar.xz files
        required=True,
        help="List of xz-compressed TAR files to process"
    )

    parser.add_argument(
        "-n", "--names",
        required=True,
        help="Text file containing list of selected FASTA file names (without extensions)"
    )

    return parser.parse_args()

# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Parse command-line arguments
    args = parse_arguments()

    selected_rnames = read_list_from_file(args.names)

    # Call the function with the list of files and selected sequence names
    extract_sequences_from_files(args.files, selected_rnames)