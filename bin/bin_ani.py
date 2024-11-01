#!/usr/bin/env python3

import pandas as pd
import argparse
import logging

def setup_logging(log_file):
    """Sets up logging to a specified log file."""
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

def bin_similarities(df, bins, bin_labels):
    """bins the ani values in the DataFrame."""
    df['ref_ani_bin'] = pd.cut(df['ani'], bins=bins, labels=bin_labels)

    bin_counts = df['ref_ani_bin'].value_counts().reindex(bin_labels)
    for bin_label, count in bin_counts.items():
        logging.info(f"Bin '{bin_label}' has {int(count)} entries.")
            
    return df

def read_tsv(file_path):
    """Reads a TSV file without headers and assigns default column names."""
    df = pd.read_csv(file_path, sep='\t', header=None, names=['query', 'reference', 'ani'])
    return df

def sample_from_bins(df, n):
    """Randomly samples n Query entries from each bin."""
    sampled_df = df.groupby('ref_ani_bin').apply(lambda x: x.sample(n=min(len(x), n)))
    sampled_df.reset_index(drop=True, inplace=True)
    return sampled_df

def save_to_tsv(df, output_path):
    """Saves the DataFrame to a TSV file."""
    df.to_csv(output_path, sep='\t', index=False)

def split_bin_list(input_string):
    """Parses a comma-separated string into a list of floats."""
    return list(map(float, input_string.split(',')))

def generate_labels(bins):
    """Generates percentage-based range labels from bins."""
    labels = []
    for i in range(len(bins) - 1):
        start = round((1 - bins[i]) * 100, 1)
        end = round((1 - bins[i + 1]) * 100, 1)
        labels.append(f"{start}-{end}%")

    return labels

def main(input_tsv, output_tsv, bin_string, n, allow_outsiders):
    bins = split_bin_list(bin_string)

    if allow_outsiders: #add a lower range to catch all
        bins = [0.0] + bins
    
    bin_labels = generate_labels(bins)
    
    df = read_tsv(input_tsv)
    
    df = bin_similarities(df, bins, bin_labels)

    # Remove rows with unassigned categories
    df = df.dropna(subset=['ref_ani_bin'])

    if n:
        df = sample_from_bins(df, n)
    
    save_to_tsv(df , output_tsv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="bin ani_similarity values from TSV")
    parser.add_argument("--input_tsv", required=True, help="Path to input TSV file of three columns: ref, query, ani")
    parser.add_argument("--output_tsv", required=True, help="Path to save output binned tsv")
    parser.add_argument("-n", type=int, help="Number of QUERY entries to sample from each bin")
    parser.add_argument("--bins", type=str, default='0.80,0.95,0.99,0.995,0.998,1', help="Comma-separated list of bin edges, e.g., '0.98,0.99,0.995,0.998,1'.")
    parser.add_argument("--assign_outsiders", action='store_true', help="add a new bin to the lower side to catch those which fall below the lowest edge")
    parser.add_argument("--log_file", type=str, default="binning.log", help="Path to save log file for binning warnings")
    
    args = parser.parse_args()

    setup_logging(args.log_file)
    
    main(args.input_tsv, args.output_tsv, args.bins, args.n, args.assign_outsiders)
