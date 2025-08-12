# Clinvar Re-Processing

## context

ClinVar is valuable resource in identifying known Pathogenic & Benign variants within a genomic dataset. By aggregating evidence from a range of submitters, we can utilise the crowd-sourced information to annotate current data with established clinical relevance.

ClinVar entries consist of:

* Individual Submissions, representing an assertion made by a submitter about the impact of an individual allele.
* Allele summaries, which produce a top-line decision about each allele by aggregating all relevant submissions. Any conflicts between submissions at the same allele are resolved in ClinVar based on [this logic](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#agg_germline).

During work on the Talos application, we found multiple known-pathogenic variants which were not identified via ClinVar due to conflicting submissions. [An example](https://ncbi.nlm.nih.gov/clinvar/variation/10/): despite 24 Pathogenic submissions to only 2 Benign, the variant is given an overall status of `Conflicting interpretations`. Whilst this is accurate based on ClinVar's internal logic, it obfuscates the bias towards pathogenicity present in the individual submissions. When we annotate a dataset with ClinVar consequences, all we have is this top-line decision, meaning that we are unable to flag such variants for more manual scrutiny.

The role of Talos is not to make clinical decisions, but to identify variants of interest for further review by analysts. In this setting we want to flag variants where manual review of the submissions could signal a variant is worth consideration, even if it doesn't appear pathogenic within the strict aggregation logic of ClinVar. To this end we created a manual re-curation of ClinVar which:

* Allows for specific submitters to be removed from consideration (i.e. so that when we run benchmarking analysis on cohorts, we can blind our ClinVar annotations to entries originating from that cohort)
* Defers to submissions after mainstream acceptance of ACMG criteria (estimated start 2016), aiming to prioritise variants which have been reviewed in light of the latest clinical guidelines
* Performs a more decisive summary, preferring a decision towards Pathogenic/Benign instead of defaulting any disagreements as `conflicting`

## process

The re-summary is rapid, and can be repeated at regular intervals, taking the latest available clinvar submissions each time it runs. The files used are the `submission_summary`and `variant_summary` present on [the NCBI clinvar FTP site](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/). We bin submissions into few discrete categories, a subset of those offered in the original ClinVar data: Benign, Conflicting, Pathogenic, Uncertain/VUS, and Unknown.

1. Iterate over all individual submissions, removing any from blacklisted providers. Collect all retained submissions per-allele.
2. For each allele, if any retained submissions were last edited after 2015 (representative ACMG date), reduce submissions to only those. If no subs are from after 2015, retain all.
3. Find a summary 'rating' across all alleles, checking these scenarios until a match is found:

   * If an Expert Review/Clinical Guideline submission is present - choose that rating.
   * If both Pathogenic and Benign submissions are present, check for a confident majority (default values: >= 60% in majority, <= 20% in minority). If there is a clear majority, choose as the overall rating.
   * If both Pathogenic and Benign subs are present, but no clear majority, assign `Conflicting`.
   * If >60% of submissions at the allele are `Uncertain`, rate as `Uncertain`.
   * If any Pathogenic submissions, take `Pathogenic`
   * If any Benign submissions, take `Benign`
   * No satisfied conditions - `Unknown/VUS`

4. A similar approach is followed for determining a `star` rating:

   * If any submissions are `Practice Guideline` -> `4 stars`
   * If any submissions are `Expert Review` -> `3 stars`
   * If any submissions `Criteria Provided` -> `1 stars`
   * Default -> `0 Stars`

At this stage we have each allele with a summary and star rating. The allele ID is matched up with the corresponding variant coordinates and ref/alt alleles from the variant summary file, then the whole object is written in multiple forms: as a Hail Table, indexed on Locus and Alleles, ready to be used in annotation within Hail; as a JSON file containing one dictionary per line, each representing one new entry; as a VCF containing all Pathogenic-rated SNVs (for use in the subsequent PM5 ACMG criteria step).

* Summarise

  * Retrieve the Submission and Variant files from NCBI's ClinVar FTP server. This is not done in code, but a bash script is provided to automate retrieval.
  * Re-summarise all submissions, saving the results as a TSV file & a Hail Table
  * Filter the Table to all Pathogenic SNVs, and export the result as a VCF

| locus        | alleles    | allele_id | clinical_significance | gold_stars |
|--------------|------------|-----------|-----------------------|------------|
| "chr1:12345" | ["A", "G"] | 789       | "Pathogenic"          | 1          |

* Annotate
  * Annotate the pathogenic SNV VCF with VEP, exporting as a TSV

* PM5 Re-Index
  * Iterate over all variants in the VCF, identifying all missense variants and their affected codon
  * Collect a lookup of the ClinVar entries relevant to each [transcript + residue] combination
  * Export this codon-indexed data, so that each [transcript + residue] combination links to all ClinVar pathogenic Missense variants affecting that codon. This is useful for PM5 ACMG criteria, which is an association of _this_ variant with _other_ missense changes at the same residue.
  * This final column has an arbitrary number of elements, consisting of a `+`-delimited String, each entry being an `AlleleID::GoldStars` string, where `AlleleID` is the unique identifier for the ClinVar allele, and `GoldStars` is the number of stars assigned to that allele.

| Transcript | Codon | ClinVar Alleles                          |
|------------|-------|------------------------------------------|
| ENST12345  | 123   | AlleleID::GoldStars                      |
| ENST12345  | 678   | AlleleID::GoldStars+AlleleID::GoldStars  |
