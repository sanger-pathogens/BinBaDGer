#  ATB dataset generator




## Usage of scripts

The ATB dataset generator relies on a number of scripts to help with the generation of appropriately filtered data and assembly clustering/deduplication.

### filter_metadata.py

This script facilitates filtering an input TSV file using a number of per-column filters (while ensuring that columns are correctly interpreted with a given type).

The column, filter and datatype are specified using the path to an input TSV manifest using the `-f` option. Valid filters are any string that can be provided to [`pandas.DataFrame.query()`](https://pandas.pydata.org/pandas-docs/version/2.2/reference/api/pandas.DataFrame.query.html) and datatypes are `int`, `float`, `datetime`, `bool` and `str` (use other types at your peril). Example manifest:

```
column	filter	datatype
center_name	center_name.str.contains("Wellcome Sanger Institute", na=False)	str
read_count	"read_count > 2500000"	int
collection_date	"2012 < collection_date"	datetime
```

By default, it will automatically remove any rows where a column value cannot be converted to the given type. For stricter validation, errors can be raised upon an invaid value (`-error_on_invalid_type` option). Missing values will be inferred, but additional values to be considered missing can be supplied with `--missing_values [MISSING_VALUES ...]`.

For the output TSV, you can select which columns you would like to output with the `--select` option, e.g. `--select sample_accession collection_date`. The header can optionally be removed from output using `--remove_header`.

Example command:
```
filter_metadata.py -f filter_manifest.tsv -i metadata.tsv -m missing not_provided not_available -s sample_accession -o test_output.tsv -r
```

You can find more help on any options using the `-h` option:
```
usage: filter_metadata.py [-h] [--input INPUT] [--filter_manifest FILTER_MANIFEST] [--select SELECT [SELECT ...]] [--missing_values MISSING_VALUES [MISSING_VALUES ...]]
                          [--remove_header] [--error_on_invalid_type] [--output OUTPUT] [--logfile LOGFILE]

Filter rows of a TSV file based on conditions.

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Path to the input TSV file. (default: None)
  --filter_manifest FILTER_MANIFEST, -f FILTER_MANIFEST
                        Path to a TSV manifest specifying the columns in the input TSV to which filter and datatype conversion should be applied. The manifest should contain 3
                        columns: column, filter, datatype. The column should match a column in the input TSV. Entries in 'column' column should be unique. The filter should match
                        a string that could be supplied to `pd.DataFrame.query()`, e.g. 'age > 30'. The datatype can be int, float, datetime, bool and str. (default: None)
  --select SELECT [SELECT ...], -s SELECT [SELECT ...]
                        Specify columns to select in the output DataFrame. By default, all columns will be selected. (default: None)
  --missing_values MISSING_VALUES [MISSING_VALUES ...], -m MISSING_VALUES [MISSING_VALUES ...]
                        Specify values that should be interpreted as missing values. (default: None)
  --remove_header, -r   Remove the header from the output TSV. (default: False)
  --error_on_invalid_type, -e
                        During type conversion, upon encountering a value in the column that cannot be converted to the given datatype, raise an error. Default behaviour is to
                        remove the row that contains an invalid value. (default: False)
  --output OUTPUT, -o OUTPUT
                        Path to the output file to save the filtered DataFrame. (default: None)
  --logfile LOGFILE, -l LOGFILE
                        Path to the log file. (default: filter_metadata-2024-49-30_10-49-53.log)
```

## Information on using sketchlib
To use these indices you will need to install the sketchlib software, which is
available from https://github.com/bacpop/sketchlib.rust. You must have the
[rust toolchain](https://www.rust-lang.org/tools/install) installed, clone the
repository, then run `cargo install --path .`. Use `sketchlib -h` to see the
help or e.g. `sketchlib dist -h` for explanation of subcommands.

The distributed index is sketch size 1024 with k=17.
You can create different ranges using the sketch subcommand.

To calculate distances of a subset of the data, use a command such as:
```
sketchlib dist -v -k 17 --subset Haemophilus_influenzae.txt --ani atb_sketchlib_v020 --threads 4 > dists.txt
```
Where the --subset file contains the list of samples you want to include.
Removing --ani will calculate Jaccard distances.

To query a set of assemblies against the index, first sketch the query sequences:
sketchlib sketch -v -o query -k 17 -f queries.tsv -s 1000
(where queries.tsv contains query samples with name and file location, tab separated)

Then query the index:
```
sketchlib dist -v -k 17 atb_sketchlib_v020 query
```

When a queryset is sketched using the sketchlib command it will return two files
.skd and .skm

Both of these are required to calculate distances down the line.