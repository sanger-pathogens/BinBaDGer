#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def load_data(input_file):
    # Load the data from TSV
    df = pd.read_csv(input_file, sep='\t', header=None, names=['Sample', 'Reference', 'ANI'])
    df['ANI'] = (df['ANI'] * 100).round(2)

    return df

def plot_histogram(df, row_count, min_ani):
    plt.figure(figsize=(10, 6))
    sns.histplot(df['ANI'], bins=20, kde=True)
    plt.title('Distribution of ANI Values')
    plt.xlabel('ANI')
    plt.ylabel('Frequency')
    plt.xlim(min_ani, 100)
    plt.figtext(0.15, 0.85, f'Number of Sample values: {row_count}', fontsize=12, ha='left')
    plt.savefig('ani_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_boxplot(df, row_count):
    plt.figure(figsize=(10, 6))
    sns.boxplot(y=df['ANI'])
    plt.title('Boxplot of ANI Values')
    plt.ylabel('ANI')
    plt.figtext(0.15, 0.85, f'Number of Sample values: {row_count}', fontsize=12, ha='left')
    plt.savefig('ani_boxplot.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_violinplot(df, row_count):
    plt.figure(figsize=(12, 8))
    sns.violinplot(y='ANI', data=df, inner='box')
    plt.title('Distribution of ANI Values')
    plt.ylabel('ANI')
    plt.figtext(0.15, 0.85, f'Number of Sample values: {row_count}', fontsize=12, ha='left')
    plt.savefig('ani_violinplot.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_heatmap(df, row_count):

    df['Sample'] = df.groupby('Sample').cumcount().astype(str).radd(df['Sample'] + '_') #little sus this needs to be added

    #This one was getting cut off

    max_height = min(50, max(8, 0.3 * len(df)))

    plt.figure(figsize=(12, max_height))
    data_map = df.pivot(index='Sample', columns='Reference', values='ANI')
    
    # Sort the data_map by the highest ANI value in each row (sample)
    data_map = data_map.loc[data_map.max(axis=1).sort_values(ascending=False).index]

    sns.heatmap(data_map)
    plt.title('Distribution of ANI Values')

    # Place figtext above the plot
    plt.figtext(0.15, 0.92, f'Number of Sample values: {row_count}', fontsize=12, ha='left')
    
    plt.savefig('ani_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

def main(input_file):
    df = load_data(input_file)
    
    row_count = df.shape[0]
    
    min_ani = df['ANI'].min()

    # Load ggplot style
    sns.set_theme(style="whitegrid")

    #plot all
    plot_histogram(df, row_count, min_ani)
    plot_boxplot(df, row_count)
    plot_violinplot(df, row_count)
    plot_heatmap(df, row_count)

if __name__ == '__main__':
    # Argument parser for command-line options
    parser = argparse.ArgumentParser(description='Visualize ANI values from a TSV file.')
    parser.add_argument('input_file', type=str, help='Input TSV file with ANI data')
    
    args = parser.parse_args()
    main(args.input_file)
