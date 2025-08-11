# ClinVar, re-summarised

## Motivation

During the creation of [Talos](https://www.github.com/populationgenomics/automated-interpretation-pipeline), a tool for identifying clinically relevant variants in large cohorts, we use [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) ratings as a contributing factor in determining pathogenicity. During development of this tool we determined that the default summaries generated in ClinVar were highly conservative; see the [table here](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#agg_germline) describing the aggregate classification logic.

## Content

This repository contains an alternative algorithm ([described here](docs/algorithm.md)) for re-aggregating the individual ClinVar submissions, generating decisions which favour clear assignment of pathogenic/benign ratings instead of defaulting to 'conflicting'. These ratings are not intended as a replacement of ClinVar's own decisions, but may provide value by showing that that though conflicting submissions exist, there is a clear bias towards either benign or pathogenic ratings.

We aim to re-run this process monthly, and publish the resulting files on Zenodo You can download this pre-generated bundle here: https://zenodo.org/records/16792026

## Primary Outputs

* Hail Table and TSV of all revised decisions
* Hail Table and TSV of all Pathogenic missense changes, indexed on Transcript and Codon. This is usable as a PM5 annotation resource.

## Usage

### Download Results

We aim to generate data monthly, and publish the results on Zenodo. The latest version of the data can be found at:

> https://zenodo.org/records/16777475

### Local Running

#### Downloading input files

A NextFlow workflow is provided to run the ClinvArbitration process locally. To use this process you will need reference files:

- a reference genome, in FASTA format
- a GFF3 file, containing gene annotations for the reference genome
- the files containing raw ClinVar submissions and variant details

A directory ([data](data)) and a script ([download_data.sh](data/download_files.sh)) are provided to download and store the required files. Running this script from the `data` directory will download and unpack all required files. The location these files are downloaded to matches the expected location in the Nextflow config, so you can run the workflow immediately after downloading.

The ClinVar Variant and Submission summary files are updated weekly. You should delete your local copy and re-download each time you run this workflow, to ensure you're capturing the latest data.

#### Running the workflow

The ClinvArbitration workflow can be run containerised, or locally. By default, the reference data will be read from a directory called `data`, and the outputs written to a directory `nextflow_outputs`.

Local execution requires:

- a Nextflow installation, to operate the workflow
- a Python environment, with the ClinvArbitration package and its dependencies installed
  - this can be actioned with `pip install .` from the root of this repository
- BCFtools, to annotate the ClinVar variants with gene information

```bash
nextflow -c nextflow/nextflow.config \
    run nextflow/clinvarbitration.nf
```

A containerised execution requires:

- a Nextflow installation, to operate the workflow
- a Docker installation, to run the workflow in a container

Step 1: build the Docker image:

```bash
docker build -t clinvarbitration:local .
```

Step 2: run the workflow using the Docker image:`

```bash
nextflow -c nextflow/nextflow.config \
    run nextflow/clinvarbitration.nf \
    -with-docker clinvarbitration:local
```

## CPG-Flow

Internally at CPG, this workflow is run using [CPG-Flow](https://github.com/populationgenomics/cpg-flow), an in-house Hail Batch based workflow executor. The following elements relate to that workflow:

* an [example config file](src/clinvarbitration/config_template.toml), with enough entries populated that a standard CPG user could dry-run the workflow locally
* a [workflow runner script](src/clinvarbitration/run_workflow.py)
* a definition of all [workflow stages](src/clinvarbitration/stages.py)

The intention is that once the Dockerfile within this repository is used, this workflow can be triggered like so:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/clinvarbitration:PR_24 \
    --config new_clinvarbitration.toml \
    --dataset seqr \
    --description 'resummarise_clinvar' \
    -o resummarise_clinvar \
    --access-level test \
    run_workflow
```

A config file is required containing a few entries, some relating to this workflow specifically, some relating to cpg-flow setup:

* `workflow.driver_image`: populated by analysis-runner, points to _this_ docker image
* `site_blacklist`: list of ClinVar submitters to ignore. Useful in removing noise, or blinding to _self_ submissions
* `ref_fasta`: required to run bcftools csq. Must match the `genome_build`
* `genome_build`: used to decide whether ClinVar/Annotation is sourced using GRCh37 or GRCh38 (default)

## Acknowledgements

* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar), for providing the data which this process is based on
