# ClinVar, reimagined

As part of work on the AIP (Automated Interpretation Pipeline) project, we are leveraging ClinVar data as one of the many factors involved in the determination of whether any individual variant is likely to be pathogenic. During development on that project, we identified a number of variants in ClinVar which had an overwhelming number of pathogenic evidence submissions overruled by a minority of benign, or uncertain. For our purposes, we are not aiming to automatically solve cases, instead surfacing a small number of variants as candidates for manual review.

We decided that a logical approach would be to take the raw data from ClinVar, and reprocess it to generate a new dataset which would be more useful for our purposes. This repository contains the code used to generate that dataset.

At CPG we leverage Hail, a python-based analysis framework which leverages Apache Spark to perform distributed computation. This repository contains a number of scripts which are designed to be run using Hail, and can either be executed once hail is installed, or using a public Hail Docker image, sourced from [DockerHub](https://hub.docker.com/r/hailgenetics/hail/tags).
