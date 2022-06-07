import ancestor_visualisation as av
import pandas as pd

def get_col_vals(df, col, val, return_col='Entry'):
    if not return_col:
        return_col = col
    cdf = df.loc[df[col] == val]
    return cdf[return_col].tolist()


anc_df = pd.read_csv(snakemake.input.csv)

tree_path = snakemake.input.tree



tree = av.load_tree(snakemake.input.tree, snakemake.input.aln)

seq_pos = (73, 83)




# Highlight the tree nodes
# tree, ts = ancestral_trace(tree, seqs_to_trace_to, seq_pos)

if snakemake.wildcards.anc_col in anc_df.columns:

    highlight_nodes = get_col_vals(anc_df, snakemake.wildcards.anc_col, snakemake.wildcards.anc_val)

else:
    print ("We couldn't find this column - {snakemake.wildcards.anc_col} in the subset")
    highlight_nodes = ['N0']




tree, ts = av.highlight_tree_nodes(tree, highlight_nodes, seq_pos)

# print (tree)

# Save the tree
print ("\nTree image has been written to " + snakemake.output.img)
tree.render(snakemake.output.img, tree_style=ts, dpi=50)