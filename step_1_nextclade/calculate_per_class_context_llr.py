"""
Consolidated per-context LLR pipeline for molnupiravir-like sequence detection.

Pipeline stages:
  1. Load nextclade output with base (type-level) LLR already computed
  2. Compute trinucleotide context for each private substitution
  3. Split into MOV-like / Normal groups using base LLR threshold (llr = 6) 
  4. Build per-context probability tables (Molnupiravir vs Normal) for each
     of the four substitution types with enough signal: G>A, A>G, C>T, T>C
  5. Apply those tables to compute a context-conditional LLR per sequence
     per substitution type
  6. Sum base LLR + all four context LLRs -> sum_llrs
  7. Threshold sum_llrs to produce final MOV call set
"""

import numpy as np
import pandas as pd
import tqdm
from collections import Counter
from scipy.stats import multinomial
from Bio import Entrez, SeqIO

# ----------------------------------------------------------------------------
# Config -- update paths before running
# ----------------------------------------------------------------------------

RAW_NEXTCLADE_TSV = "/Users/reem/Mov/final_results_2026.tsv"  
CLASS_DISTRIBUTION_TSV = "/Users/reem/Downloads/estimated_mutation_distribution.tsv"  # pre-known 12-class MutationType distribution, used for mutational class LLR
OUTPUT_DIR = "/Users/reem/Mov/nextclade_results/"

SUBSTITUTION_TYPES = ["G>A", "A>G", "C>T", "T>C"]  # the four types modeled at context resolution
LLR_THRESHOLD = 6    # Empirical relaxed threshold to avoid missing any MOV-like sequences in BTE analysi later on.
BASES = ["A", "C", "G", "T"]


# ----------------------------------------------------------------------------
# Stage 1: load data with base (type-level) LLR already present
# ----------------------------------------------------------------------------

def load_raw_data(path):
    """Load raw nextclade output."""
    iterator = pd.read_csv(
        path,
        sep="\t",
        usecols=["seqName", "privateNucMutations.unlabeledSubstitutions"],
        dtype={"seqName": str, "privateNucMutations.unlabeledSubstitutions": str},
        chunksize=1000,
    )
    return pd.concat([chunk for chunk in tqdm.tqdm(iterator, desc="Loading raw data")])


def compute_class_llr(df, class_distribution_path=CLASS_DISTRIBUTION_TSV):
    """
    Mutational class LLR using the expected 12-class mutation distribution under molnupiravir vs normal models.
    pM/pN here represent MutationType only (no context yet).
    """
    probs_df = pd.read_csv(class_distribution_path, sep="\t")
    mutation_types = probs_df["MutationType"].str.replace("→", ">").tolist()
    pM = probs_df["Molnupiravir"].to_numpy(dtype=float)
    pN = probs_df["Normal"].to_numpy(dtype=float)
    assert np.isclose(pM.sum(), 1, atol=1e-2), f"Molnupiravir class probs sum to {pM.sum():.4f}, not ~1"  # Make sure the probabilities are valid
    assert np.isclose(pN.sum(), 1, atol=1e-2), f"Normal class probs sum to {pN.sum():.4f}, not ~1"

    def llr_for_row(counts_dict):
        counts = np.array([counts_dict.get(mt, 0) for mt in mutation_types], dtype=float)
        n = counts.sum()
        if n == 0:
            return np.nan
        llM = multinomial.logpmf(counts, n=n, p=pM)
        llN = multinomial.logpmf(counts, n=n, p=pN)
        return float(llM - llN)

    df["LLR"] = df["Counts"].apply(llr_for_row)
    return df


# ----------------------------------------------------------------------------
# Stage 2: reference genome + context annotation
# ----------------------------------------------------------------------------

def fetch_reference_genome(accession="NC_045512.2", email="theo@theo.io"):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)


def get_mut_type(mut_string):
    """e.g. 'A234G' -> 'A>G'"""
    if isinstance(mut_string, str) and len(mut_string) >= 2:
        return mut_string[0] + ">" + mut_string[-1]
    return ""


def get_context(genome_seq, mutation):
    """Trinucleotide context of a private substitution, e.g. 'A234G' -> 'xAy'"""
    if isinstance(mutation, str) and len(mutation) >= 2:
        pos = int(mutation[1:-1]) - 1  # nextclade positions are 1-based
        return genome_seq[pos - 1:pos + 2]
    return ""


def build_spectrum(subs, contexts):
    """Combine mutation type + context, e.g. 'C[G>A]T'"""
    if not isinstance(subs, str) or not isinstance(contexts, str):
        return ""
    if not subs.strip() or not contexts.strip():
        return ""
    subs_list = subs.split(",")
    ctx_list = contexts.split(",")
    out = []
    for mutation, context in zip(subs_list, ctx_list):
        if len(context) >= 2 and mutation:
            out.append(f"{context[0]}[{mutation}]{context[-1]}")
    return ",".join(out)


def annotate_contexts(df, reference_genome):
    """Convenience wrapper if you want Stage 2 context annotation standalone
    (main() inlines these same steps so compute_class_llr can run in between)."""
    df["subs"] = df["privateNucMutations.unlabeledSubstitutions"].apply(
        lambda x: ",".join(get_mut_type(item) for item in x.split(",")) if isinstance(x, str) else ""
    )
    df["Counts"] = df["subs"].apply(
        lambda x: dict(Counter(x.split(","))) if isinstance(x, str) and x.strip() else {}
    )
    df["context"] = df["privateNucMutations.unlabeledSubstitutions"].apply(
        lambda x: ",".join(get_context(reference_genome, item) for item in x.split(",")) if isinstance(x, str) else ""
    )
    df["spectrum"] = df.apply(lambda row: build_spectrum(row["subs"], row["context"]), axis=1)
    return df


def generate_all_possible_contexts(sub):
    """All 16 trinucleotide contexts for a given substitution type, e.g. 'G>A' -> ['A[G>A]A', ...]"""
    return [f"{a}[{sub}]{b}" for a in BASES for b in BASES]


def count_contexts_for_type(spectrum_str, sub_type):
    """Count contexts for one specific substitution type within a sequence's full spectrum string."""
    counts = Counter()
    if not isinstance(spectrum_str, str):
        return {}
    for mut in spectrum_str.split(","):
        if mut[2:5] == sub_type:
            counts[mut] += 1
    return dict(counts)


def annotate_type_context_counts(df, substitution_types=SUBSTITUTION_TYPES):
    for sub in substitution_types:
        col = sub.replace(">", "to") + "_counts"  # e.g. GtoA_counts -- ">" avoided per VERIFY #3
        df[col] = df["spectrum"].apply(lambda s, sub=sub: count_contexts_for_type(s, sub))
    return df


# ----------------------------------------------------------------------------
# Stage 3: build per-context probability tables from labeled groups
# ----------------------------------------------------------------------------

def get_mean_context_proportions(df, llr_condition, counts_col, llr_col="LLR"):
    """Mean per-context proportion of a substitution type, within sequences passing llr_condition."""
    df_filtered = df[llr_condition(df[llr_col])].copy()

    def to_proportions(counts_dict):
        total = sum(counts_dict.values())
        if total == 0:
            return {}
        return {k: v / total for k, v in counts_dict.items()}

    props = df_filtered[counts_col].apply(to_proportions)
    pivot = pd.json_normalize(props).fillna(0)
    pivot = pivot.loc[(pivot != 0).any(axis=1)]  # drop sequences with zero of this substitution type
    return pivot.mean()


def build_prob_tables_once(df, output_dir=OUTPUT_DIR, substitution_types=SUBSTITUTION_TYPES,
                            group_threshold=LLR_THRESHOLD, rebuild=False):
    """
    Builds and saves the per-context probability tables only if they don't
    already exist on disk (or if rebuild=True is passed explicitly). This
    enforces "build once, freeze forever" -- confirmed as the intended
    behavior. Returns the loaded tables.
    """
    tables = {}
    for sub in substitution_types:
        fname = f"{output_dir}{sub.replace('>', 'to')}probs_new.tsv"
        if not rebuild and pd.io.common.file_exists(fname):
            print(f"Loading existing frozen table: {fname}")
            tables[sub] = pd.read_csv(fname, sep="\t")
        else:
            print(f"Building new table (rebuild={rebuild}, exists={pd.io.common.file_exists(fname)}): {fname}")
            table = build_prob_table(df, sub, group_threshold=group_threshold)
            table.to_csv(fname, sep="\t", index=False)
            tables[sub] = table
    return tables


def build_prob_table(df, sub_type, llr_col="LLR", group_threshold=LLR_THRESHOLD):
    """
    Builds Molnupiravir vs Normal mean-proportion table for one substitution
    type's 16 contexts.
    """
    counts_col = sub_type.replace(">", "to") + "_counts"
    all_contexts = generate_all_possible_contexts(sub_type)

    mov_probs = get_mean_context_proportions(df, lambda x: x > group_threshold, counts_col, llr_col)
    normal_probs = get_mean_context_proportions(df, lambda x: x <= group_threshold, counts_col, llr_col)

    mov_probs = mov_probs.reindex(all_contexts, fill_value=0)
    normal_probs = normal_probs.reindex(all_contexts, fill_value=0)

    df_prob = pd.DataFrame({
        "Mutational_Context": all_contexts,
        "Molnupiravir": mov_probs.values,
        "Normal": normal_probs.values,
    })

    #atol is absolute tolerance, set yto 1e-2 to allow for small numerical errors in floating point summation 
    mov_sum = df_prob["Molnupiravir"].sum()
    normal_sum = df_prob["Normal"].sum()
    assert np.isclose(mov_sum, 1, atol=1e-2) or mov_sum == 0, \
        f"[{sub_type}] Molnupiravir proportions sum to {mov_sum:.4f}, not ~1 -- decomposition is invalid, do not proceed"
    assert np.isclose(normal_sum, 1, atol=1e-2) or normal_sum == 0, \
        f"[{sub_type}] Normal proportions sum to {normal_sum:.4f}, not ~1 -- decomposition is invalid, do not proceed"

    return df_prob


# ----------------------------------------------------------------------------
# Stage 4: apply probability tables to compute per-context LLR
# ----------------------------------------------------------------------------

def calculate_llr(count_dict, pM, pN, contexts):
    counts = np.array([count_dict.get(ctx, 0) for ctx in contexts], dtype=float)
    n = counts.sum()
    if n == 0:
        return np.nan
    llM = multinomial.logpmf(counts, n=n, p=pM)
    llN = multinomial.logpmf(counts, n=n, p=pN)
    return float(llM - llN)


def apply_context_llrs(df, prob_tables, substitution_types=SUBSTITUTION_TYPES):
    """prob_tables: dict of sub_type -> df_prob (from build_prob_table or loaded from disk)"""
    llr_cols = []
    for sub in substitution_types:
        counts_col = sub.replace(">", "to") + "_counts"
        llr_col = sub.replace(">", "to") + "_llr"
        llr_cols.append(llr_col)
        df_prob = prob_tables[sub]
        contexts = df_prob["Mutational_Context"].values.tolist()
        pM = df_prob["Molnupiravir"].values.tolist()
        pN = df_prob["Normal"].values.tolist()
        df[llr_col] = df[counts_col].apply(lambda cd: calculate_llr(cd, pM, pN, contexts))
    return df, llr_cols


# ----------------------------------------------------------------------------
# Stage 5: combine and call
# ----------------------------------------------------------------------------

def combine_and_call(df, llr_cols, base_llr_col="LLR", threshold=LLR_THRESHOLD):
    # NaN -> 0 assumes "no evidence from this substitution type" contributes nothing,
    # rather than being missing data (i.e. that substitution type just didn't
    # occur in this sequence's private mutations). This matches the confirmed
    # design: base_llr_col ("LLR") is always the primary/base signal, the four
    # context LLRs are additional refinement on the MOV-relevant types only.
    df["sum_llrs"] = df[base_llr_col].fillna(0) + sum(df[c].fillna(0) for c in llr_cols)
    df_mov = df[df["sum_llrs"] > threshold].copy()
    return df, df_mov

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main(rebuild_tables=False):
    # Stage 1: raw data + class LLR from the pre-known
    # class distribution -- this is the "LLR" column, computed here, not loaded.
    df = load_raw_data(RAW_NEXTCLADE_TSV)
    df["subs"] = df["privateNucMutations.unlabeledSubstitutions"].apply(
        lambda x: ",".join(get_mut_type(item) for item in x.split(",")) if isinstance(x, str) else ""
    )
    df["Counts"] = df["subs"].apply(
        lambda x: dict(Counter(x.split(","))) if isinstance(x, str) and x.strip() else {}
    )
    df = compute_class_llr(df)

    # Stage 2: trinucleotide context annotation (needed for the 4 context tables)
    reference_genome = fetch_reference_genome()
    df["context"] = df["privateNucMutations.unlabeledSubstitutions"].apply(
        lambda x: ",".join(get_context(reference_genome, item) for item in x.split(",")) if isinstance(x, str) else ""
    )
    df["spectrum"] = df.apply(lambda row: build_spectrum(row["subs"], row["context"]), axis=1)
    df = annotate_type_context_counts(df)

    # Stage 3: load probability tables (build only if missing / rebuild_tables=True).
    prob_tables = build_prob_tables_once(df, rebuild=rebuild_tables)

    # Stage 4/5: apply tables, combine, call
    df, llr_cols = apply_context_llrs(df, prob_tables)
    df, df_mov = combine_and_call(df, llr_cols)

    print(f"Sequences called MOV-like (sum_llrs > {LLR_THRESHOLD}): {len(df_mov)} / {len(df)}")

    drop_cols = ["privateNucMutations.unlabeledSubstitutions", "subs", "Counts",
                 "context", "spectrum"] + [sub.replace(">", "to") + "_counts" for sub in SUBSTITUTION_TYPES]
    df_mov.drop(columns=[c for c in drop_cols if c in df_mov.columns]).to_csv(
        f"{OUTPUT_DIR}mov_2026.tsv", sep="\t", index=False
    )

    keep_cols = ["seqName", "LLR"] + llr_cols + ["sum_llrs"]
    df[keep_cols].to_csv(f"{OUTPUT_DIR}final_llrs_2026.tsv", sep="\t", index=False)

    print("Done.")


if __name__ == "__main__":
    main()