






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