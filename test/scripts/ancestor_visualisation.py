from Bio.Seq import Seq
from Bio import SeqIO
from ete3 import (
    PhyloTree,
    TreeStyle,
    TextFace,
    add_face_to_node,
    SeqMotifFace,
    SequenceFace,
    NodeStyle,
    faces,
    ImgFace,
    CircleFace,
    AttrFace,
)


def load_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = PhyloTree(tree_path, alignment=aln_path, format=1, alg_format="fasta")
    return tree


def processable_node(node, highlight_nodes):
    """
    If there is still a node in highlight nodes that is a descendant of current node, we want to keep going
    """
    check = False
    for name in highlight_nodes:
        if name in [x.name for x in node.get_descendants()]:
            return True
        else:
            check = False
    return check


def highlight_tree_nodes(
    tree,
    highlight_nodes=["N0"],
    node_colour="red",
    node_style="sphere",
    node_size=10,
    leaf_colour="green",
    leaf_style="sphere",
    leaf_size=10,
    seq_pos=False,
    seq_motif=False,
    display_all_nodes=False,
    seq_change=False,
):
    ts = TreeStyle()
    # disable default PhyloTree Layout
    ts.layout_fn = lambda x: True

    for n in tree.traverse():
        # If the node is either in highlight nodes or one of its descendents is
        if (
            processable_node(n, highlight_nodes)
            or n.name in highlight_nodes
            or display_all_nodes
        ):
            if n.name in highlight_nodes or display_all_nodes:
                col_name = n.name

                if not n.is_leaf():
                    N = AttrFace("name", fsize=24, fgcolor="black")
                    n.add_face(N, 1, position="branch-top")

                    nstyle = NodeStyle()
                    nstyle["shape"] = node_style
                    nstyle["fgcolor"] = node_colour
                    nstyle["size"] = node_size
                    nstyle["hz_line_type"] = 2
                    n.set_style(nstyle)

                elif n.is_leaf():
                    nstyle = NodeStyle()
                    nstyle["shape"] = leaf_style
                    nstyle["fgcolor"] = leaf_colour
                    nstyle["size"] = leaf_size
                    nstyle["hz_line_type"] = 2
                    n.set_style(nstyle)

                if seq_pos:
                    subseq = [n.sequence[x] for x in seq_pos]

                    parentseq = (
                        [n.up.sequence[x] for x in seq_pos] if not n.is_root() else ""
                    )

                    # If we're just tracking where the sequence changes, don't display the same subsequence
                    if seq_change and "".join(subseq) == "".join(parentseq):
                        pass
                    else:
                        if seq_motif:
                            # Get all the subseqs under this point
                            childseqs = [[subseq]]
                            for child in n.get_descendants():
                                print(child)
                                childseqs.append([child.sequence[x] for x in seq_pos])
                            #                             childseqs = [[child.sequence[x] for child in n.get_descendants()] for x in seq_pos]

                            S = SeqMotifFace(childseqs, seq_format="seq", width=20)
                            n.add_face(S, 1, position="branch-bottom")

                        else:
                            S = SequenceFace(
                                subseq, interactive=True, col_w=20, fsize=18
                            )
                            n.add_face(S, 1, position="branch-bottom")

        else:
            # We've reached a node which isn't one of the nodes to highlight and doesn't have a node to highlight as one of its descendents

            if n.name not in highlight_nodes:
                format_text = " extant sequence" if len(n) == 1 else " extant sequences"
                N = TextFace(" " + str(len(n)) + format_text, fsize=24, fgcolor="black")
                n.add_face(N, 1, position="branch-right")
                wid = max(len(n) * 0.3, 20)
                n.add_face(ImgFace("./triangle.png", width=wid), 0)
                n.img_style["draw_descendants"] = False
                n.name = ""

    return tree, ts


def ancestral_trace(tree, seqs_to_trace_to, seq_pos, seq_change=False):
    """
    Given a tree, extant sequences, and a position, perform an ancestal trace
    """
    highlight_nodes = set()
    for seq_name in seqs_to_trace_to:
        node = tree & seq_name
        ancestor_list = node.get_ancestors()
        highlight_nodes.update(x.name for x in ancestor_list)
        highlight_nodes.add(seq_name)

    tree, ts = highlight_tree_nodes(
        tree, highlight_nodes=highlight_nodes, seq_change=seq_change, seq_pos=seq_pos
    )

    return tree, ts
