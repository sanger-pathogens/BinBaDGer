#!/usr/bin/env python3

import pandas as pd
import argparse

def bin_similarities(df, bins, bin_labels):
    """bins the ani values in the DataFrame."""
    df['ref_ani_bin'] = pd.cut(df['ani'], bins=bins, labels=bin_labels)
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

def main(input_tsv, output_tsv, bin_string, n):
    bins = split_bin_list(bin_string)
    
    bin_labels = generate_labels(bins)
    
    df = read_tsv(input_tsv)
    
    df = bin_similarities(df, bins, bin_labels)

    if n:
        df = sample_from_bins(df, n)
    
    save_to_tsv(df , output_tsv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="bin ani_similarity values from TSV")
    parser.add_argument("--input_tsv", required=True, help="Path to input TSV file of three columns: ref, query, ani")
    parser.add_argument("--output_tsv", required=True, help="Path to save output binned tsv")
    parser.add_argument("-n", type=int, help="Number of QUERY entries to sample from each bin")
    parser.add_argument("--bins", type=str, default='0.80,0.95,0.99,0.995,0.998,1', help="Comma-separated list of bin edges, e.g., '0.98,0.99,0.995,0.998,1'.")
    
    args = parser.parse_args()
    
    main(args.input_tsv, args.output_tsv, args.bins, args.n)
