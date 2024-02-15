#!/usr/bin/env bash

set -ex

# create a docker image from this repository
docker build -t hail_clinvar:example --platform linux/amd64 .

# make local copies of the NCBI data files required as input using wget
# create a directory called data, if one doesn't already exist
if [ ! -d data ]; then
    mkdir data
fi
wget -O data/variant_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
wget -O data/submission_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz

# run the docker image to generate the summarised output
docker run --platform linux/amd64 -v "$(pwd)/data":/data hail_clinvar:example \
    /bin/bash -c "python3 /scripts/resummarise.py -v /data/variant_summary.txt.gz -s /data/submission_summary.txt.gz -o /data/clinvar_summary"

# upon completion, this will have generated files in the data directory:
# - data/clinvar_summary.json - a JSON file containing the summarised data entries, one json object per line
# - data/clinvar_summary_for_VEP.vcf.bgz - a bgzipped VCF which can be used in VEP annotation
# - data/clinvar_summary.vcf.bgz - a bgzipped file containing the pathogenic SNV entries in VCF format
# - data/clinvar_summary.ht - a Hail Table containing the summarised data entries

# This is where you should run VEP on data/clinvar_summary.vcf.bgz, with protein consequence annotation per transcript
# Let's imagine you did that, and the result is in data/pathogenic_annotated.vcf.bgz
# I've enclosed a 10-variant example of this, as annotated by https://www.ensembl.org/Homo_sapiens/Tools/VEP
docker run --platform linux/amd64 -v "$(pwd)/data":/data hail_clinvar:example \
    /bin/bash -c "python3 /scripts/clinvar_by_codon_from_vcf.py -i /data/pathogenic_annotated.vcf.bgz -o /data/pm5"

# upon completion, this will generate files in the data directory:
# - data/pm5.json - a JSON file containing the PM5 results, one JSON object per line
# - data/pm5.ht - a Hail table containing the same PM5 results
