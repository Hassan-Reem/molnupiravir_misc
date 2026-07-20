# Step 1: Nextclade Processing & LLR Calculation

This step runs Nextclade on raw sequences from GISAID, calculates per-class and per-context
log-likelihood ratios (LLRs) for molnupiravir-like mutational signatures,
and visualizes the resulting G→A / C→T mutational spectra in comparison to non-molnupiravir sequences.

## Pipeline order

1. `nextclade_runner.py` — runs Nextclade on input sequences (`batch.fa`)
2. `calculate_per_class_context_llr.py` — parses Nextclade output, computes
   trinucleotide context, calculates base-class + per-context LLRs, and
   applies the probability tables below to produce final LLR calls
3. `plot_gtoa_ctot_spectrum.ipynb` — plots G→A and C→T mutational spectra
   for MOV-like vs Normal groups.

## Probability tables

- `GtoA_probs.tsv`, `AtoG_probs.tsv`, `CtoT_probs.tsv`, `TtoC_probs.tsv` —
  per-context (16 trinucleotide contexts) Molnupiravir vs Normal
  probability tables, built once from the original labeled MOV/Normal
  split on a threshold of LLR = 6 and used for all subsequent LLR calculations.
  
  **Do not regenerate these** they are the mean per-context proportions, computed separately for two groups of sequences:
  Molnupiravir column: sequences with class LLR ≥ 6 ("MOV-like").
  Normal column: sequences with base LLR < 6 ("Normal").
  The entire downstream analysis (and the
  published results) rely on these probability tables.


## Outputs
- `final_llrs_26.tsv` — per-sequence LLR values (used downstream in step_2_nextclade and Usher_analysis)

## Exploratory / archived
See `exploratory/` for earlier drafts and one-off analyses not used in the
final pipeline.