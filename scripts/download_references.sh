#!/usr/bin/env bash

module load ncbi-datasets-cli/14.27.0

accessions_file="${1}"  # File of RefSeq (GCF*) or GenBank (GCA*) assembly accessions
outdir="${2}"  # Output directory where you want to root the directory tree of references

accessions=()
while read accession; do
  accessions+=("${accession}")
done < "${accessions_file}"


datasets download genome accession "${accessions[@]}" --include genome,seq-report
tmp_dir=$(mktemp -d)
unzip -o -d "${tmp_dir}" ncbi_dataset.zip

readarray -t dl_accessions < <(jq -r '.accession' ${tmp_dir}/ncbi_dataset/data/assembly_data_report.jsonl)
readarray -t sanitized_species < <(jq -r '.organism.organismName' ${tmp_dir}/ncbi_dataset/data/assembly_data_report.jsonl | cut -d ' ' -f 1,2 | tr " " "_" | tr [:upper:] [:lower:])

for i in "${!dl_accessions[@]}"; do
  accession="${dl_accessions[i]}"
  species="${sanitized_species[i]}"
  accession_dir="${outdir}/${species}/${accession}"
  mkdir -p "${accession_dir}"
  mv ${tmp_dir}/ncbi_dataset/data/${accession}/*.{fna,jsonl} ${accession_dir}
done

rm -r ncbi_dataset.zip "${tmp_dir}"

