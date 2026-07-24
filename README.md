# Molnupiravir-Associated Mutagenesis in SARS-CoV-2

Analysis pipeline for detecting and characterising molnupiravir-associated mutational signatures in global SARS-CoV-2 genomic data using a probabilistic log-likelihood ratio (LLR) framework.

## Overview

Molnupiravir is a mutagenic antiviral that induces a characteristic pattern of transition mutations in SARS-CoV-2. This repository contains the code and analysis workflows used to identify molnupiravir-associated mutational signatures across large-scale SARS-CoV-2 genomic datasets and to investigate their evolutionary and epidemiological consequences.

The pipeline combines:

- Nextclade-based sequence analysis
- Mutation-class and trinucleotide-context likelihood modelling
- A multinomial log-likelihood ratio (LLR) framework
- UShER phylogenetic analysis
- Big Tree Explorer (BTE) for large-scale phylogenetic traversal
- Temporal and geographic inference using MetaFitch
- Benchmarking against previously published molnupiravir-detection methods
- Analysis of transmission dynamics, fitness, recurrent mutations and geographic distribution

## Repository Structure

### `step_1_nextclade/`

Initial sequence processing and mutational signature modelling.

This step includes:

- Running Nextclade on batches of SARS-CoV-2 sequences
- Extracting private mutations
- Characterising the molnupiravir-associated mutation spectrum
- Calculating per-mutation-class and trinucleotide-context likelihood ratios
- Generating the empirical probability tables used by the LLR framework

### `step_2_nextclade/`

Threshold selection and identification of sequences requiring addition to the UShER phylogeny.

This step includes:

- Screening Nextclade outputs using empirical LLR thresholds
- Identifying sequences absent from the UShER tree
- Selecting candidate sequences for subsequent addition to the tree using the UShER placement pipeline.

### `UShER_analysis/`

Large-scale phylogenetic analysis of molnupiravir-associated lineages.

This step includes:

- Traversing the UShER mutation-annotated phylogeny using Big Tree Explorer (BTE)
- Applying the LLR framework to phylogenetic nodes
- Inferring temporal and geographic metadata using MetaFitch
- Benchmarking model performance using Australia and France as positive and negative epidemiological reference settings, respectively
- Selecting classification thresholds using mutation-class and trinucleotide-context LLRs
- Comparing the LLR framework with previously published threshold-based methods
- Analysing geographic and temporal distributions
- Assessing transmission and cluster size
- Evaluating fitness consequences
- Identifying recurrent mutations
- Investigating the 2025 United States molnupiravir-associated cluster
- Notebooks for generating the figures used in the manuscript.

## Analysis Workflow

The overall analysis proceeds through the following stages:

1. **Sequence processing and quality control** using Nextclade
2. **Mutational signature characterisation** and construction of empirical probability models
3. **LLR calculation** for mutation classes and trinucleotide contexts
4. **Screening of sequences absent from the UShER phylogeny**
5. **Large-scale phylogenetic classification** of molnupiravir-associated nodes
6. **Temporal and geographic inference** using MetaFitch
7. **Benchmarking and comparison** with existing detection methods
8. **Evolutionary and epidemiological analysis** of identified lineages
9. **Manuscript figure generation**

## Reproducibility

The repository contains the analysis scripts used to generate the results presented in the associated manuscript. The underlying SARS-CoV-2 genomic sequences and associated metadata are not included because their redistribution is restricted by the terms and conditions governing access to the GISAID database.

The analyses require access to the relevant GISAID datasets and external phylogenetic resources. Researchers wishing to reproduce the analysis should obtain access to the appropriate datasets directly through the relevant data providers and comply with their applicable terms of use.
Where possible, intermediate analysis outputs, derived probability tables and summary data that do not contain restricted underlying sequence data are provided in the repository.

Please refer to the individual scripts and accompanying documentation for required input files, software dependencies and analysis parameters.

## Contact

For questions regarding the analysis or repository, please contact:

**Reem Hassan**  
London School of Hygiene & Tropical Medicine  
reem.hassan@lshtm.ac.uk
