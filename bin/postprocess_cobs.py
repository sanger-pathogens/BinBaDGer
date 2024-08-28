#!/usr/bin/env python3

import argparse
import os
import re

def get_nb_kmers(cobs_line):
    p = cobs_line.split("\t")
    nb_kmers = int(p[-1])
    return nb_kmers

def remove_rnd_id(cobs_line):
    _, _, r = cobs_line.partition("_")
    return "_" + r

def process_cobs_output(hits_to_keep, files):
    identifier = ""
    total_count = 0
    combined_data = []
    
    # Process each file
    for file in files:
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    # Handle the header
                    header_parts = line.strip().split("\t")
                    file_identifier = header_parts[0]
                    count = int(header_parts[1])
                    
                    # Set the identifier if not already set
                    if not identifier:
                        identifier = file_identifier
                    
                    # Accumulate the count
                    total_count += count
                else:
                    # Process the remaining lines
                    combined_data.append(remove_rnd_id(line))

    # Write the new header
    print(f"{identifier}\t{total_count}")

    # Initialize tracking variables
    min_kmers = 0
    i = 0

    # Sort combined data by kmers count (assuming tab-separated and last column is the count)
    sorted_data = sorted(combined_data, key=get_nb_kmers, reverse=True)
    
    # Output the sorted data while keeping only the top N hits
    for y in sorted_data:
        i += 1
        if i <= hits_to_keep:
            print(y, end="")
            if i == hits_to_keep:
                min_kmers = get_nb_kmers(y)
        elif get_nb_kmers(y) == min_kmers:
            print(y, end="")

def main():
    parser = argparse.ArgumentParser(
        description="Postprocess cobs output: keep top n hits (+ties), remove random identifiers, and combine files"
    )

    parser.add_argument(
        '-n',
        metavar='int',
        dest='keep',
        required=True,
        type=int,
        help='Number of best hits to keep',
    )

    parser.add_argument(
        'files',
        metavar='F',
        type=str,
        nargs='+',
        help='List of files to process',
    )

    args = parser.parse_args()

    process_cobs_output(args.keep, args.files)

if __name__ == "__main__":
    main()
