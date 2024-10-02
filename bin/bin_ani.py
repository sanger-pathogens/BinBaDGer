#!/usr/bin/env python3

import pandas as pd
import argparse

def bin_similarities(df, bins, bin_labels):
    """Bins the ANI values in the DataFrame."""
    df['Bin'] = pd.cut(df['ANI'], bins=bins, labels=bin_labels)
    return df

def read_tsv(file_path):
    """Reads a TSV file without headers and assigns default column names."""
    df = pd.read_csv(file_path, sep='\t', header=None, names=['Query', 'Reference', 'ANI'])
    return df

def save_to_tsv(df, output_path):
    """Saves the DataFrame to a TSV file."""
    df.to_csv(output_path, sep='\t', index=False)

def sample_from_bins(df, n):
    """Randomly samples n Query entries from each bin."""
    sampled_df = df.groupby('Bin').apply(lambda x: x.sample(n=min(len(x), n)))
    sampled_df.reset_index(drop=True, inplace=True)
    return sampled_df

def main(input_tsv, output_tsv, n):
    bins = [0.98, 0.99, 0.995, 0.998, 1]
    bin_labels = ['2%', '1%', '0.5%', '0.2%']
    
    df = read_tsv(input_tsv)
    
    df = bin_similarities(df, bins, bin_labels)

    sampled_df = sample_from_bins(df, n)
    
    save_to_tsv(sampled_df , output_tsv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bin ANI_SIMILARITY values from TSV")
    parser.add_argument("input_tsv", help="Path to input TSV file")
    parser.add_argument("output_tsv", help="Path to save output TSV file")
    parser.add_argument("n", type=int, help="Number of QUERY entries to sample from each bin")
    
    args = parser.parse_args()
    
    main(args.input_tsv, args.output_tsv, args.n)
