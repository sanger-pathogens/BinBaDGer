#!/usr/bin/env python3

import argparse
import os
import re

def get_kmer_count(cobs_line):
    columns = cobs_line.split("\t")
    kmer_count = int(columns[-1])
    return kmer_count

def remove_random_id(cobs_line):
    _, _, remainder = cobs_line.partition("_")
    return remainder

def process_cobs_output(hits_to_keep, input_files, show_count):
    main_identifier = ""
    total_hit_count = 0
    combined_lines = []
    identifiers_list = []
    
    # Process each file
    for file_path in input_files:
        with open(file_path, 'r') as file:
            for line_index, line in enumerate(file):
                if line_index == 0:
                    # Handle the header
                    header_parts = line.strip().split("\t")
                    file_identifier = header_parts[0]
                    hit_count = int(header_parts[1])
                    
                    # Set the main identifier if not already set
                    if not main_identifier:
                        main_identifier = file_identifier
                    
                    # Accumulate the hit count
                    total_hit_count += hit_count
                    
                    # Store the identifier and optionally the count
                    if show_count:
                        identifiers_list.append(f"{file_identifier}\t{hit_count}")
                    else:
                        identifiers_list.append(file_identifier)
                else:
                    # Process the remaining lines
                    combined_lines.append(remove_random_id(line))

    # Initialize tracking variables
    min_kmer_count = 0
    line_counter = 0

    # Sort combined lines by kmer count (assuming tab-separated and last column is the count)
    sorted_lines = sorted(combined_lines, key=get_kmer_count, reverse=True)
    
    # Output the sorted data while keeping only the top N hits
    for sorted_line in sorted_lines:
        line_counter += 1
        if line_counter <= hits_to_keep:
            if show_count:
                print(sorted_line, end="")
            else:
                # Remove the kmer count before printing
                line_without_count = "\t".join(sorted_line.split("\t")[:-1])
                print(line_without_count)
            if line_counter == hits_to_keep:
                min_kmer_count = get_kmer_count(sorted_line)
        elif get_kmer_count(sorted_line) == min_kmer_count:
            if show_count:
                print(sorted_line, end="")
            else:
                line_without_count = "\t".join(sorted_line.split("\t")[:-1])
                print(line_without_count)


def main():
    parser = argparse.ArgumentParser(
        description="Postprocess COBS output: keep top N hits (+ties), remove random identifiers, and combine files"
    )

    parser.add_argument(
        '-n',
        metavar='int',
        dest='hits_to_keep',
        required=True,
        type=int,
        help='Number of top hits to keep',
    )

    parser.add_argument(
        '--show_count',
        action='store_true',
        help='Show the hit count for each identifier',
    )

    parser.add_argument(
        'input_files',
        metavar='F',
        type=str,
        nargs='+',
        help='List of input files to process',
    )

    args = parser.parse_args()

    process_cobs_output(args.hits_to_keep, args.input_files, args.show_count)

if __name__ == "__main__":
    main()
