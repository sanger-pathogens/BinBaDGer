# Running BinBaDGer on the Sanger farm

The AllTheBacteria COBS indexes, Sketchlib database and assemblies are stored in lustre and set as defaults for the `--cobs_base`, `--sketchlib_db` and `--assembly_base` parameters respectively. There are also Kraken2 databases pre-downloaded from [here](https://benlangmead.github.io/aws-indexes/k2) in `/data/pam/software/kraken2` that you can use with the `--kraken_db` parameter- the default is `/data/pam/software/kraken2/k2_standard_16gb_20240904`.

You will need to load the Nextflow and Singularity modules:

```
module load nextflow
module load ISG/singularity
```

Then you can run the pipeline with bsub. We recommend using the oversubscribed queue and requesting 4GB memory, for example:

```
bsub -q oversubscribed -J binbadger -R "select[mem>4000] rusage[mem=4000]" -M4000 -o binbadger.%J.o -e binbadger.%J.e nextflow run path/to/repo/main.nf --manifest manifest.csv
```