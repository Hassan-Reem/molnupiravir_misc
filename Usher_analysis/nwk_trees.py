import toytree
import toyplot
import os


"""This script visualizes Newick phylogenetic trees, highlighting MOV samples (matches.txt) as red nodes,
   and saves the resulting trees as SVG files for easy viewing."""


# Load Newick files
mov_samples = "/Users/reem/matches.txt"
newick_dir = "/Users/reem/nwk_files"
output_dir = "/Users/reem/nwk_trees"

with open(mov_samples) as f:
    mov_samples = [line.strip() for line in f]


for nwk_file in os.listdir(newick_dir):
    filename = os.fsdecode(nwk_file)
    if filename.endswith(".nw"):
        newick_path = os.path.join(newick_dir, filename)
        tree = toytree.tree(newick_path)

    node_colors = ["blue"] * tree.nnodes

    # Color leaves: blue for default, red if in matches
    for leaf in tree[:tree.ntips]:  # slice of all leaf nodes
        if leaf.name in mov_samples:
            node_colors[leaf.idx] = "red"
        else:
            node_colors[leaf.idx] = "blue"
    canvas, axes, marks = tree.draw(
        node_sizes=13, 
        node_mask=False, 
        node_colors=node_colors,
        scale_bar=True,
        width=800,
        height=800,
        edge_widths=2)


    canvas.style.update({"background-color": "white"})
    toytree.save(canvas, os.path.join(output_dir, f"{filename}.svg"))
    print(f"Saved {filename}")




    