import scripts.annot_functions as an
import scripts.ancestor_visualisation as av
import pandas as pd
import scripts.seqcurate as sc
from Bio import AlignIO


print("Adding custom annotations from ancestors")


align_df = pd.read_csv(snakemake.input.csv)

anc_df = sc.get_sequence_df(snakemake.input.aln, alignment=True, ancestor=True)

aln = AlignIO.read(snakemake.input.aln, format="fasta")

aln_dict = {seq.name: str(seq.seq) for seq in aln}

align_df.reset_index(inplace=True, drop=True)
anc_df.reset_index(inplace=True, drop=True)

frames = [align_df, anc_df]
merge_df = pd.concat(frames)

merge_df.to_csv(snakemake.output.csv, index=False)

# # Check if the KARI EC number is in the dataset we're currently looking at
# if "ec_1_1_1_86" in snakemake.wildcards.dataset:
#
#     print ("Adding KARI specific annotations")
#
#
#     # tree_path = "./GRASP_ancestors.nwk"
#     # aln_path = "./90_1e-75_motif_confirmed_long_N_C_removed_missing_motif_removed_truncated_DSE_extension_removed_further_edits_cropped_w_ancestors_cleaned.aln"
#     nodes = "N0 N1"
#
#     out_path = "./KARI_example_output.png"
#     out_path2 = "./KARI_example_output2.png"
#
#     # Load tree
#     tree = av.load_tree(snakemake.input.tree, snakemake.input.aln)
#
#     print (tree)
#     # Get the nodes of interest
#     node_list = [x for x in nodes.split()]
#
#     seq_pos = [25, 26, 27]
#
#     # Highlight the tree nodes
#     tree, ts = av.highlight_tree_nodes(tree, node_list, seq_pos=seq_pos, seq_motif=True, display_all_nodes=True)
#
#     # Save the tree
#     # tree.render(out_path2, tree_style=ts)
#
#     tree.render(out_path2, tree_style=ts, w = int(200), h = int(300), dpi=int(6000))
#     tree.render(out_path, tree_style=ts, dpi=int(6000))

# Create an ancestral trace image
