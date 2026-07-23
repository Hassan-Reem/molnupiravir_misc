"""
BTE (Big Tree Explorer) node-level LLR calculation and MOV cluster detection.

This is the tree-based counterpart to calculate_per_class_context_llr.py (step_1):
instead of per-sequence LLRs, this computes LLR per internal tree NODE (i.e.
per branch/mutation set), then identifies contiguous "MOV-like" clusters in
the UShER tree topology.

Pipeline stages:
  1. Load the annotated UShER tree (.pb) and extract per-node mutation data
  2. Compute base (class) LLR and per-context LLRs per node
  3. [EXTERNAL STEP -- see run_metafitch_externally() below] run metafitch
     for ancestral date/country reconstruction, OUTSIDE this script
  4. Merge metafitch output back onto the node table
  5. Compute descendant leaf counts per node
  6. Identify MOV cluster roots and cluster sizes -- run as a SENSITIVITY
     SWEEP across sum_llrs thresholds
  7. Save outputs, plot cluster size distributions per threshold
"""

import re
from collections import Counter

import bte
import numpy as np
import pandas as pd
from scipy.stats import multinomial
from Bio import Entrez, SeqIO

# ----------------------------------------------------------------------------
# Config -- update paths before running
# ----------------------------------------------------------------------------

TREE_PB = "/Users/reem/2026_updated_tree.pb"
TREE_NEWICK = "/Users/reem/2026_tree.nwk"  # generated from TREE_PB via matUtils extract, required by metafitch
CLASS_DISTRIBUTION_TSV = "/Users/reem/Downloads/estimated_mutation_distribution.tsv"
CONTEXT_PROBS_DIR = "/Users/reem/Mov/nextclade_results/molnupiravir_misc/step_1_nextclade/" 
METAFITCH_DATE_TSV = "/Users/reem/metafitch_date_output.tsv"
METAFITCH_COUNTRY_TSV = "/Users/reem/metafitch_country_output.tsv"
OUTPUT_DIR = "/Users/reem/demo/"

SUBSTITUTION_TYPES = ["G>A", "A>G", "C>T", "T>C"]


# ----------------------------------------------------------------------------
# Stage 1: load tree, extract per-node mutation data
# ----------------------------------------------------------------------------

def load_tree(tree_pb_path):
    return bte.MATree(tree_pb_path)
 
 
def iter_node_batches(tree, batch_size=100000):
    batch = []
    for node in tree.depth_first_expansion():
        muts = tree.get_node(node.id).mutations
        if len(muts) == 0:
            continue
        batch.append({"node_id": node.id, "mutations": muts})
        if len(batch) >= batch_size:
            yield pd.DataFrame(batch)
            batch = []
    if batch:
        yield pd.DataFrame(batch)
 
 
def fetch_reference_genome(accession="NC_045512.2", email="theo@theo.io"):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)
 
 
# ----------------------------------------------------------------------------
# Stage 2: LLR calculation (base class + per-context), reusing step_1's
# probability tables. Mutation-string parsing differs slightly from step_1
# because BTE stores mutations as a list-like column, not a nextclade string.
# ----------------------------------------------------------------------------
 
def clean_mutation_token(mut_string):
    """Strip stray quote characters left over from list-to-string conversion."""
    if not isinstance(mut_string, str):
        return ""
    return mut_string.strip().strip("'\"")
 
 
def get_mut_type(mut_string):
    mut_string = clean_mutation_token(mut_string)
    if len(mut_string) >= 3:
        return f"{mut_string[0]}>{mut_string[-1]}"
    return ""
 
 
def get_context(genome_seq, mutation):
    mutation = clean_mutation_token(mutation)
    if len(mutation) < 3 or mutation[0] not in "ACGT" or mutation[-1] not in "ACGT":
        return ""
    middle = mutation[1:-1]
    if not middle.isdigit():
        return ""
    pos = int(middle) - 1
    return genome_seq[pos - 1:pos + 2]
 
 
def build_spectrum(subs, contexts):
    if not isinstance(subs, str) or not isinstance(contexts, str):
        return ""
    if not subs.strip() or not contexts.strip():
        return ""
    subs_list = [s.strip() for s in subs.split(",")]
    ctx_list = [c.strip() for c in contexts.split(",")]
    out = []
    for mutation, context in zip(subs_list, ctx_list):
        if len(context) >= 2 and mutation:
            out.append(f"{context[0]}[{mutation}]{context[-1]}")
    return ",".join(out)
 
 
def count_contexts_for_type(spectrum_str, sub_type):
    counts = Counter()
    if not isinstance(spectrum_str, str):
        return {}
    for mut in spectrum_str.split(","):
        if mut[2:5] == sub_type:
            counts[mut] += 1
    return dict(counts)
 
 
def calculate_llr(count_dict, pM, pN, categories):
    counts = np.array([count_dict.get(c, 0) for c in categories], dtype=float)
    n = counts.sum()
    if n == 0:
        return np.nan
    llM = multinomial.logpmf(counts, n=n, p=pM)
    llN = multinomial.logpmf(counts, n=n, p=pN)
    return float(llM - llN)
 
 
def load_llr_reference_tables(class_distribution_path=CLASS_DISTRIBUTION_TSV,
                               context_probs_dir=CONTEXT_PROBS_DIR,
                               substitution_types=SUBSTITUTION_TYPES):
    probs_df = pd.read_csv(class_distribution_path, sep="\t")
    mut_types = probs_df["MutationType"].str.replace("→", ">").tolist()
    pM_class = probs_df["Molnupiravir"].to_numpy(dtype=float)
    pN_class = probs_df["Normal"].to_numpy(dtype=float)
 
    context_tables = {}
    for sub in substitution_types:
        table_path = f"{context_probs_dir}{sub.replace('>', 'to')}_probs.tsv"  # VERIFY: matches step_1's renamed convention
        prob_table = pd.read_csv(table_path, sep="\t")
        categories = prob_table["Mutational_Context"].tolist()
        pM_ctx = prob_table["Molnupiravir"].to_numpy(dtype=float)
        pN_ctx = prob_table["Normal"].to_numpy(dtype=float)
        context_tables[sub] = (pM_ctx, pN_ctx, categories)
 
    return mut_types, pM_class, pN_class, context_tables
 
 
def annotate_node_llrs_batch(df_batch, reference_genome, mut_types, pM_class, pN_class,
                              context_tables, substitution_types=SUBSTITUTION_TYPES):
    df_batch = df_batch.copy()
    df_batch["mutations_str"] = df_batch["mutations"].astype(str).str.strip("[]").str.replace("'", "", regex=False)
    df_batch["subs"] = df_batch["mutations_str"].apply(
        lambda x: ",".join(mut for mut in (get_mut_type(item) for item in x.split(",")) if mut) if isinstance(x, str) else ""
    )
    df_batch["Counts"] = df_batch["subs"].apply(
        lambda x: dict(Counter(x.split(","))) if isinstance(x, str) and x.strip() else {}
    )
    df_batch["LLR"] = df_batch["Counts"].apply(lambda cd: calculate_llr(cd, pM_class, pN_class, mut_types))
 
    df_batch["context"] = df_batch["mutations_str"].apply(
        lambda x: ",".join(mut for mut in (get_context(reference_genome, item) for item in x.split(",")) if mut) if isinstance(x, str) else ""
    )
    df_batch["spectrum"] = df_batch.apply(lambda row: build_spectrum(row["subs"], row["context"]), axis=1)
 
    llr_cols = []
    for sub in substitution_types:
        counts_col = sub.replace(">", "to") + "_counts"
        llr_col = sub.replace(">", "to") + "_llr"
        llr_cols.append(llr_col)
 
        counts_series = df_batch["spectrum"].apply(lambda s, sub=sub: count_contexts_for_type(s, sub))
        pM_ctx, pN_ctx, categories = context_tables[sub]
        df_batch[llr_col] = counts_series.apply(lambda cd, pM=pM_ctx, pN=pN_ctx, cats=categories: calculate_llr(cd, pM, pN, cats))
        del counts_series  # dict column, discard immediately once the LLR is computed from it
 
    df_batch["sum_contexts"] = df_batch[llr_cols].fillna(0).sum(axis=1)
    df_batch["sum_llrs"] = df_batch["LLR"].fillna(0) + df_batch["sum_contexts"]
 
    keep_cols = ["node_id", "LLR"] + llr_cols + ["sum_contexts", "sum_llrs"]
    return df_batch[keep_cols].copy()
 
 
# ----------------------------------------------------------------------------
# Stage 3: EXTERNAL STEP -- metafitch (ancestral date/country reconstruction)
# ----------------------------------------------------------------------------
 
def run_metafitch_externally():
    """
    THIS FUNCTION DOES NOT RUN METAFITCH. It documents the external,
    manual terminal steps that must be run before continuing this pipeline.
    metafitch performs ancestral state reconstruction (Fitch parsimony) on
    the UShER tree to infer date/country for internal nodes that don't have
    this metadata directly (only tip/leaf sequences have observed
    dates/countries; internal nodes need it inferred from their
    descendants). Run TWICE -- once for date, once for country -- producing
    two separate output files that get merged together below.
 
    Step 1 -- convert the protobuf tree to Newick (metafitch requires .nwk,
    not .pb):
 
        matUtils extract -i 2026_updated_tree.pb -t 2026_tree.nwk
 
    Step 2 -- run metafitch once per trait:
 
        # Date run
        python3 /Users/reem/miniconda3/envs/bte/lib/python3.12/site-packages/metafitch/metafitch.py \\
            -t 2026_tree.nwk \\
            -m metafitch_dates.tsv \\
            -o ./metafitch_date_output.tsv \\
            -a
 
        # Country run
        # Note: -a (allow ambiguous dates) is a date-specific flag; likely a
        # no-op here since this run has no date field, but included for
        # consistency with the date run above. Confirm whether metafitch
        # ignores it silently or whether it should be dropped for this run.
        python3 /Users/reem/miniconda3/envs/bte/lib/python3.12/site-packages/metafitch/metafitch.py \\
            -t 2026_tree.nwk \\
            -m metafitch_countries.tsv \\
            -o ./metafitch_country_output.tsv \\
            -a
 
    Inputs:
      -t  2026_tree.nwk              -- Newick tree from step 1
      -m  metafitch_dates.tsv /
          metafitch_countries.tsv    -- metadata file to insert (leaf-level
                                         observed date/country per strain,
                                         joined on "strain")
      -a                              -- allow ambiguous dates (e.g. "2022",
                                         "2023" rather than a full date) to
                                         be used in the reconstruction,
                                         instead of discarding them
 
    Outputs: METAFITCH_DATE_TSV, METAFITCH_COUNTRY_TSV (see config at top of
    this script). Both are merged together (on "strain") and then onto the
    node-level LLR table in merge_metafitch_output() below.
 
    Requires the `bte` conda environment specifically
    (/Users/reem/miniconda3/envs/bte/) -- this is a separate environment
    from wherever calculate_llrs_consolidated.py is normally run, worth
    noting explicitly in the repo README so a future run doesn't fail with
    a confusing ModuleNotFoundError the same way the Homebrew/conda
    interpreter mismatch did earlier.
    """
    raise NotImplementedError(
        "matUtils + metafitch must be run manually in a terminal (in the "
        "'bte' conda environment) -- see this function's docstring for the "
        "exact commands. This function is a placeholder marking where in "
        "the pipeline that step belongs."
    )
 
 
def merge_metafitch_output(df_nodes, date_path=METAFITCH_DATE_TSV, country_path=METAFITCH_COUNTRY_TSV):
    """Merges the two separate metafitch runs (date, country) into one
    table, then merges that onto the node-level LLR table.
 
    Confirmed: both metafitch outputs use "strain" as the join-key column
    name (renamed to "node_id" below to match df_nodes).
    """
    date_df = pd.read_csv(date_path, sep="\t")
    country_df = pd.read_csv(country_path, sep="\t")
 
    metafitch_out = pd.merge(date_df, country_df, on="strain", how="outer", suffixes=("_date_run", "_country_run"))
    metafitch_out = metafitch_out.rename(columns={"strain": "node_id"})
 
    merged = pd.merge(df_nodes, metafitch_out, on="node_id", how="left")
    return merged
 
 
# ----------------------------------------------------------------------------
# Stage 4: tree topology lookups -- parent/descendant counts, etc. These are 
# used in the MOV cluster detection below.
# ----------------------------------------------------------------------------
 
def get_parent_id(node):
    """
    Returns the parent's node_id as a string, regardless of whether
    node.parent is a MATNode object, some other node wrapper, or already a
    bare string ID.
    """
    parent = node.parent
    if parent is None:
        return None
    return parent.id if hasattr(parent, "id") else parent
 
 
def get_descendant_leaves(root_node, tree):
    leaves = []
    stack = [root_node]
    while stack:
        node = stack.pop()
        if len(node.children) == 0:
            leaves.append(node.id)
            continue
        for child in node.children:
            child_node = tree.get_node(child.id) if hasattr(child, "id") else tree.get_node(child)
            if child_node is not None:
                stack.append(child_node)
    return leaves
 
 
def get_descendant_count(node_id, tree):
    node = tree.get_node(node_id)
    if node is None:
        return None
    return len(get_descendant_leaves(node, tree))
 
 
def clean_mutations_list(x):
    if pd.isna(x):
        return []
    x = str(x).strip().strip("[]")
    return [tok.strip().strip("'\"") for tok in x.split(",") if tok.strip().strip("'\"")]
 
 
# ----------------------------------------------------------------------------
# Stage 5: identify MOV clusters using LLR and sum_contexts thresholds
# ----------------------------------------------------------------------------

LLR_THRESHOLD = 3
CONTEXT_THRESHOLD = 2 

def identify_mov_clusters(merged, tree):
    """
    A cluster root is a node whose sum_llrs exceeds `threshold` but whose
    parent does not -- this avoids double-counting nested MOV-like nodes
    within the same cluster.
    """
    mov_ids = set(merged.loc[(merged["LLR"] > LLR_THRESHOLD) & (merged["sum_contexts"] > CONTEXT_THRESHOLD), "node_id"])
    print(f"Identified {len(mov_ids)} MOV-like nodes exceeding thresholds (LLR>{LLR_THRESHOLD}, sum_contexts>{CONTEXT_THRESHOLD})")

    mov_cluster_root_ids = []
    for node_id in mov_ids:
        node = tree.get_node(node_id)
        if node is None:
            continue
        parent_id = get_parent_id(node)
        if parent_id is None or parent_id not in mov_ids:
            mov_cluster_root_ids.append(node_id)
    print(f"MOV cluster roots detected: {len(mov_cluster_root_ids):,}")

    cluster_sizes = []
    for node_id in mov_cluster_root_ids:
        node = tree.get_node(node_id)
        if node is None:
            cluster_sizes.append(0)
            continue
        n_leaves = len(get_descendant_leaves(node, tree))
        cluster_sizes.append(n_leaves
                             )
    mov_cluster_roots_df = pd.DataFrame({
        "node_id": mov_cluster_root_ids,
        "cluster_size": cluster_sizes,
    })
    return mov_cluster_roots_df
 
 
def plot_cluster_size_distribution(cluster_df, output_path=None):
    import matplotlib.pyplot as plt
    from collections import Counter as _Counter
 
    size_counts = _Counter(cluster_df["cluster_size"])
    sizes = list(size_counts.keys())
    counts = list(size_counts.values())
 
    plt.figure(figsize=(10, 5))
    bars = plt.bar(sizes, counts, width=1.0, edgecolor="black")
    plt.xlabel("# of descendant leaves", fontsize=12)
    plt.ylabel("# MOV clusters", fontsize=12)
    plt.title(f"MOV Cluster Size Distribution", fontsize=14)
    plt.gca().bar_label(bars, label_type="edge", color="black")
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path)
    plt.show()
 
 
# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
 
def get_mutations_for_node_ids(tree, node_ids):
    return {nid: clean_mutations_list(str(tree.get_node(nid).mutations)) for nid in node_ids}
 
 
def main(batch_size=100000, resume_from_llrs = False):
    tree = load_tree(TREE_PB)
    reference_genome = fetch_reference_genome()
    mut_types, pM_class, pN_class, context_tables = load_llr_reference_tables()
 
    llr_output_path = f"{OUTPUT_DIR}bte_nodes_llrs.tsv"

    if resume_from_llrs:
        print(f"Resuming from existing {llr_output_path} -- skipping batched LLR calculation.")
        df_nodes = pd.read_csv(llr_output_path, sep="\t")
    else:
        reference_genome = fetch_reference_genome()
        mut_types, pM_class, pN_class, context_tables = load_llr_reference_tables()
        first_batch = True
        n_processed = 0
        for batch_df in iter_node_batches(tree, batch_size=batch_size):
            result = annotate_node_llrs_batch(
                batch_df, reference_genome, mut_types, pM_class, pN_class, context_tables
            )
            result.to_csv(llr_output_path, sep="\t", mode="w" if first_batch else "a",
                            header=first_batch, index=False)
            n_processed += len(result)
            print(f"Processed {n_processed} nodes so far...")
            first_batch = False
            del batch_df, result  
 
  
    df_nodes = pd.read_csv(llr_output_path, sep="\t")
 
    # --- Stage 3 checkpoint: run metafitch manually here before continuing ---
    # See run_metafitch_externally() docstring for the command that needs
    merged = merge_metafitch_output(df_nodes)
    merged["num_descendants"] = merged["node_id"].apply(lambda nid: get_descendant_count(nid, tree))
    merged_output_path = f"{OUTPUT_DIR}merged_bte_final.tsv"
    merged.to_csv(merged_output_path, sep="\t", index=False)
    print("Saved merged table at:", merged_output_path)

    merged = pd.read_csv(f"{OUTPUT_DIR}merged_bte_final.tsv", sep="\t")
    clusters = identify_mov_clusters(merged, tree)
    clusters_output_path = f"{OUTPUT_DIR}bte_mov_clusters.tsv"
    clusters.to_csv(clusters_output_path, sep="\t", index=False)
    print("Saved MOV cluster table at:", clusters_output_path)

    plot_cluster_size_distribution(clusters, output_path=f"{OUTPUT_DIR}MOV_cluster_size_distribution.svg")
 
    print("Done.")
 
 
if __name__ == "__main__":
    main()
               