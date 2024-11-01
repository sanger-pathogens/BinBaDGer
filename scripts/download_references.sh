#!/usr/bin/env bash

module load ncbi-datasets-cli/14.27.0

accessions_file="${1}"  # File of RefSeq (GCF*) or GenBank (GCA*) assembly accessions
outdir="${2}"  # Output directory where you want to root the directory tree of references


# Read accessions into array 
readarray -t accessions < "${accessions_file}"

# Download data using NCBI datasets CLI
datasets download genome accession "${accessions[@]}" --include genome,seq-report
tmp_dir=$(mktemp -d)
unzip -o -d "${tmp_dir}" ncbi_dataset.zip

# Order of downloaded accessions not the same as those provided to datasets CLI, so read into array again from metadata
readarray -t dl_accessions < <(jq -r '.accession' ${tmp_dir}/ncbi_dataset/data/assembly_data_report.jsonl)
readarray -t sanitized_species < <(jq -r '.organism.organismName' ${tmp_dir}/ncbi_dataset/data/assembly_data_report.jsonl | cut -d ' ' -f 1,2 | tr " " "_" | tr [:upper:] [:lower:])

# Keep metadata
metadata_dir="${outdir}/metadata"
mkdir -p "${metadata_dir}"
mv ${tmp_dir}/ncbi_dataset/data/assembly_data_report.jsonl "${metadata_dir}"/assembly_data_report_$(date -I).jsonl

# Move references into dir structure
for i in "${!dl_accessions[@]}"; do
  accession="${dl_accessions[i]}"
  species="${sanitized_species[i]}"
  accession_dir="${outdir}/${species}/${accession}"
  mkdir -p "${accession_dir}"
  mv ${tmp_dir}/ncbi_dataset/data/${accession}/*.{fna,jsonl} ${accession_dir}
done

rm -r ncbi_dataset.zip "${tmp_dir}"

