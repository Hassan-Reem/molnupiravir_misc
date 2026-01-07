import re
import os
from pathlib import Path
import glob
import SeqIO
from matutils import matutils

folder = Path("/Users/reem/usher_output/")

for filepath in glob.glob(os.path.join(folder, "subtree-*.txt")):
    print(filepath)

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("node_"):
                continue

            matutils extract -i "/Users/reem/Downloads/usher_pruned.pb" -s 

