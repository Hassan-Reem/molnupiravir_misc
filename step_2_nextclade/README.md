# Step 2: Identify missing sequences & LLR Threshold Selection

Identifies GISAID FASTA sequences absent from the UShER phylogenetic tree and
prepares them for placement, and empirically selects the sum_llrs
threshold used to call MOV-like sequences for subsequent addition to trhe UShER tree.

## Files

1. `map_missing_seqs_for_usher_placement.py` extracts FASTA sequences
   not present in the UShER tree (based on `unmatched.txt`) and outputs
   them in batches for subsequent UShER placement.
2. `select_llr_threshold.ipynb` empirically selects the sum_llrs
   threshold (= 6) used for MOV-like calling, by comparing detection
   counts in Australia (MOV-associated) vs. France (non-MOV-associated)
   sequences across a range of thresholds. Threshold 6 was chosen to
   preserve sensitivity (538 Australia detections vs. 243 at t=10),
   accepting a small increase in background signal (30 France detections
   vs. 0 at t=10), since additional specificity is provided by
   downstream BTE-based filtering.

## Data (gitignored)
Intermediate `.tsv` files are excluded from
version control. Regenerate via the scripts above as needed.