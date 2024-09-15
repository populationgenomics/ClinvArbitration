# ClinVar, reimagined

## motivation

During the creation of [Talos](https://www.github.ocom/populationgenomics/automated-interpretation-pipeline), a tool for identifying clinically relevant variants in large cohorts, we leveraged the entries in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) as a contributing factor in determining pathogenicity. During development of this tool we determined that the default summaries generated in ClinVar were highly conservative; see the [table here](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#agg_germline) describing the aggregate classification logic.

## content

This repository contains an alternative algorithm ([described here](docs/algorithm.md)) for re-visiting the ClinVar default results, generating alternative decisions which tend towards clear assignment of pathogenic/benign ratings. These ratings are not intended as a replacement of ClinVar's own decisions, but may provide value in situations where an analyst would benefit from seeing that though conflicting submissions exist, there is a bias towards either benign or pathogenic ratings.

## outputs

Our intention with this repository is to make this code and process available, as well as periodically producing releases containing the resulting data files for consumption in other analyses.

This currently generates a few key outputs:

* JSON file of all revised decisions
* Hail Table of all revised decisions
* VCF of all revised decisions; this can be used as a custom annotation source in VEP
* VCF of all Pathogenic-SNVs, for annotation & feeding into the second stage; ACMG criteria PM5 analysis

[this script](example_script.sh) shows the steps involved in generating the data indicated above. The VCF form can be used as a custom annotation source in VEP. See this syntax available in VEP >= 110:

```bash
./vep [...] --custom file=clinvar_for_VEP.vcf.gz,short_name=CPG_ClinVar,format=vcf,type=exact,coords=0,fields=allele_id%gold_stars%clinical_significance
```

## acknowledgements

* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar), for providing the data which this process is based on

---

At CPG we leverage Hail, a python-based analysis framework which leverages Apache Spark to perform distributed computation. This repository contains a number of scripts which were initially designed to be run using Hail, and can either be executed once hail is installed, or using a public Hail Docker image, sourced from [DockerHub](https://hub.docker.com/r/hailgenetics/hail/tags). A Dockerfile included in this repository will build a custom image capable of locally executing all scripts.

Our aim with this repository is to carry out a periodic reprocessing of ClinVar data, and to make the results available to the wider community in a range of formats. We hope that this will be useful to others who are working with ClinVar data, and that it will be a useful resource for those who are working on similar projects.
