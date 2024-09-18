#!/usr/bin/env python3

import argparse
import random

def get_kmer_count(cobs_line):
    """Extract the k-mer count from a line."""
    columns = cobs_line.split("\t")
    return int(columns[-1])

def remove_random_id(cobs_line):
    """Remove the random identifier from a line."""
    _, _, remainder = cobs_line.partition("_")
    return remainder

def process_cobs_file(file_path):
    """Process the file to collect identifiers and lines with valid hits."""
    identifiers = set()
    lines = []

    with open(file_path, 'r') as cobs_file:
        for line in cobs_file:
            line = line.strip()
            if line.startswith('*'):
                # Collect identifiers
                header_parts = line.split("\t")
                file_identifier = header_parts[0]
                identifiers.add(file_identifier)
            else:
                # Process valid lines
                lines.append(remove_random_id(line))
    
    return identifiers, lines

def sort_and_print_lines(lines, hits_to_keep, show_count, output_file):
    """Sort lines by k-mer count and print the top N hits."""
    sorted_lines = sorted(lines, key=get_kmer_count, reverse=True)
    selected_lines = sorted_lines[:hits_to_keep]

    write_lines_to_file(selected_lines, show_count, output_file)

def random_and_print_lines(lines, hits_to_keep, show_count, output_file, seed=None):
    """Randomly select and print up to N hits."""
    
    if seed is not None:
        random.seed(seed)
    
    selected_lines = random.sample(lines, hits_to_keep)
    
    # Print the selected lines
    write_lines_to_file(selected_lines, show_count, output_file)

def step_and_print_lines(lines, hits_to_keep, show_count, output_file, seed):
    """Systematically select and print up to N hits across the entire range."""
    
    # Set the seed for reproducibility if provided
    if seed is not None:
        random.seed(seed)

    # Determine the step size and starting index
    step = total_lines // hits_to_keep
    start = random.randint(0, step - 1)
    selected_indices = [start + i * step for i in range(hits_to_keep)]

    # Select the lines based on the calculated indices
    selected_lines = [lines[idx] for idx in selected_indices]

    # Print the selected lines
    write_lines_to_file(selected_lines, show_count, output_file)

def write_lines_to_file(selected_lines, show_count, output_file):
    seen_lines = set()

    with open(output_file, 'w') as out:
        for line in selected_lines:
            if show_count:
                out.write(line)  # Writes the full line including the count
            else:
                # Remove the kmer count before writing
                line_without_count = "\t".join(line.split("\t")[:-1])
                
                if line_without_count not in seen_lines:
                    out.write(line_without_count + "\n")
                    seen_lines.add(line_without_count)

def process_cobs_output(hits_to_keep, input_files, show_count, selection_method, seed):
    """Process the COBS output files to filter, optionally sort, and select hits."""
    combined_lines = []

    output_file = 'combined_matches.txt'

    for file_path in input_files:
        # Collect identifiers and process lines from the file
        _, file_lines = process_cobs_file(file_path)
        combined_lines.extend(file_lines)
    
    total_lines = len(combined_lines)

    print(f'total found {total_lines} before subset to {hits_to_keep}')

    if total_lines <= hits_to_keep or selection_method == 'top':
        # Sort and print lines
        sort_and_print_lines(combined_lines, hits_to_keep, show_count, output_file)

    elif selection_method == 'random':
        # Randomly select and print lines
        random_and_print_lines(combined_lines, hits_to_keep, show_count, output_file, seed)

    elif selection_method == 'stepwise':
        # Step over the data and print lines across the entire range
        step_and_print_lines(combined_lines, hits_to_keep, show_count, output_file, seed)
    else:
        raise Exception(f"no method for {selection_method}")

def main():
    parser = argparse.ArgumentParser(
        description="Postprocess COBS output: keep top N hits (+ties), remove random identifiers, and combine files"
    )

    parser.add_argument(
        '-n',
        metavar='int',
        dest='hits_to_keep',
        default=500,
        type=int,
        help='Number of top hits to keep, default: 500 if set to much higher number there could not be enough possible matches, \
        also if the number of hits to keep is greater than the number of matches\
        all matches are returned with no selection',
    )

    parser.add_argument(
        '--show_count',
        action='store_true',
        help='Show the hit count for each identifier',
    )

    parser.add_argument(
        '--selection_method',
        metavar='method',
        type=str,
        choices=['top', 'random', 'stepwise'],
        default='top',
        help='Method to select the hits: "top" (default), "random", or "stepwise"',
    )

    parser.add_argument(
        '--seed',
        metavar='int',
        type=int,
        default=None,
        help='Seed for random selection (for reproducibility)',
    )

    parser.add_argument(
        'input_files',
        metavar='F',
        type=str,
        nargs='+',
        help='List of input files to process',
    )

    args = parser.parse_args()

    process_cobs_output(args.hits_to_keep, args.input_files, args.show_count, args.selection_method, args.seed)

if __name__ == "__main__":
    main()
