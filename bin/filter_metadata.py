#!/usr/bin/env python3

import argparse
import pandas as pd
import logging
import sys

from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter rows of a TSV file based on conditions.')
    
    parser.add_argument('--input', '-i', type=Path, help='Path to the input TSV file.')
    parser.add_argument('--cleaning_filters', '-f1', type=str, nargs='+', help="Filter conditions in the format 'column operator value'. Example: 'age > 30'.")
    parser.add_argument('--logic_filters', '-f2', type=str, nargs='+', help="Filter conditions in the format 'column operator value'. Example: 'age > 30'.")
    parser.add_argument('--column_dtypes', '-d', type=str, nargs='+', help="Force the datatypes of columns. Specify the column and datatype using the syntax 'col:type'.")
    parser.add_argument('--output', '-o', type=str, help='Path to the output file to save the filtered DataFrame.')
    
    return parser.parse_args()


def build_query_string(conditions):
    """ Converts filter conditions into a valid query string for pandas """
    query_list = []
    for condition in conditions:
        query_list.append(condition)
    return ' and '.join(query_list)


def clean_dataframe(df):
    """ Remove leading and trailing whitespaces from column names and values in DataFrame """
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

    return df

def apply_filters(df, filters: list[str]):
    query_string = build_query_string(filters)
    try:
        filtered_df = df.query(query_string)
    except Exception as e:
        raise e
        # logging.error(f"Error in filtering: {e}")
        # sys.exit(1)
    return filtered_df

def parse_column_type_list(column_types: list):
    """
    Parses a list of "column:type" mappings and returns a dictionary with column names as keys and types as values.
    
    Parameters:
    column_type_list (list of str): A list of column-to-type mappings in the format "column:type".
    
    Returns:
    dict: A dictionary where keys are column names and values are target types.
    """
    parsed_column_types = {}
    for column_type in column_types:
        col, dtype = column_type.split(":")
        parsed_column_types[col] = dtype
    return parsed_column_types

def safe_convert_column(df, column, dtype):
    """
    Attempts to convert a DataFrame column to the specified data type.
    Removes rows where any value prevents conversion.
    
    Parameters:
    df (pd.DataFrame): The DataFrame containing the column to convert.
    column (str): The name of the column to convert.
    dtype (str): The target data type (e.g., 'int', 'float', 'datetime', etc.).
    
    Returns:
    pd.DataFrame: The DataFrame with the column converted, and rows with invalid values removed.
    """
    if dtype == 'int' or dtype == 'float':
        # Use pd.to_numeric with errors='coerce' for numeric conversion
        df[column] = pd.to_numeric(df[column], errors='coerce')
    elif dtype == 'datetime':
        # Use pd.to_datetime with errors='coerce' for datetime conversion
        df[column] = pd.to_datetime(df[column], errors='coerce')
    else:
        # For other types, try using astype with errors handling
        try:
            df[column] = df[column].astype(dtype)
        except ValueError as e:
            logging.error(f"Error converting column '{column}' to type '{dtype}': {e}")

    # Drop rows where conversion resulted in NaN (failed conversions)
    df = df.dropna(subset=[column])

    return df

def apply_column_types(df, column_types: dict):
    """
    Applies the given column-to-type mappings to the DataFrame.
    
    Parameters:
    df (pd.DataFrame): The DataFrame to modify.
    column_types (dict): A dictionary where keys are column names and values are target types.
    
    Returns:
    pd.DataFrame: The DataFrame with updated column types.
    """
    for col, dtype in column_types.items():
        if col in df.columns:
            df = safe_convert_column(df, col, dtype)
        else:
            logging.error(f"Column '{col}' not found in DataFrame.")
    
    return df
    #     try:
    #         # Convert the column to the specified dtype
    #         df[col] = df[col].astype(dtype)
    #         print(df[col].dtype)
    #     except KeyError:
    #         print(f"Column '{col}' not found in DataFrame.")
    #     except ValueError as e:
    #         raise e
    #         # print(f"Error converting column '{col}' to {dtype}: {e}")
    
    # return df


def main():
    args = parse_arguments()
    column_types = parse_column_type_list(args.column_dtypes)
    
    
    # Read and clean the metadata TSV
    df = pd.read_csv(args.input, sep='\t')
    # df = clean_dataframe(df)
    # Ensure we still encode missing values correctly
    # df.replace('', pd.NA, inplace=True)
    
    # TODO With safe conversion of dtype we may no longer need to apply cleaning filters!
    if args.cleaning_filters:
        df = apply_filters(df, args.cleaning_filters)
    
    # Not necessary if we use a cleaning filter like "collection_data.notnull()"
    # df = df.dropna(subset=["collection_date"])

    # Found problematic values in collection date column:
    # print(df.iloc[1053]["collection_date"]) # not applicable
    # print(df.iloc[1116]["collection_date"])  # missing
    # print(df.iloc[1126]["collection_date"])  # not provided


    df = df.convert_dtypes()  # TODO Still necessary?
    df = apply_column_types(df, column_types)
    
    if args.logic_filters:
        df = apply_filters(df, args.logic_filters)


    # Output the filtered DataFrame
    if args.output:
        df.to_csv(args.output, sep='\t', index=False)
        logging.info(f"Filtered data saved to {args.output}")
    else:
        # Print to stdout
        print(df)

if __name__ == "__main__":
    main()
