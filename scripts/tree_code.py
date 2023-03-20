from ete3 import (
    PhyloTree,
    TreeStyle,
    TextFace,
    add_face_to_node,
    SeqMotifFace,
    NodeStyle,
    faces,
    ImgFace,
    CircleFace,
    AttrFace,
)
import random


def load_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = PhyloTree(tree_path, alignment=aln_path, format=1, alg_format="fasta")
    return tree


def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=5)
        faces.add_face_to_node(N, node, 0, position="aligned")


def get_example_tree(tree, color_dict=None, annot_dict=None, col=None, rotate=None):
    used_colours = set()
    for n in tree.traverse():
        if rotate and n.name in rotate:
            print("swapping")
            n.swap_children()
        if n.is_leaf():
            if n.name in annot_dict:
                n.img_style["bgcolor"] = color_dict[annot_dict[n.name]]

                used_colours.add(annot_dict[n.name])

            ts = TreeStyle()
            ts.layout_fn = layout
            ts.show_leaf_name = False
            ts.mode = "c"
            # ts.arc_start = -180  # 0 degrees = 3 o'clock
            # ts.arc_span = 180
            # ts.root_opening_factor = 1

            # print (used_colours)

            for k, v in color_dict.items():
                # Check that the specific colour was actually present in the tree we're annotating
                if k in used_colours:
                    ts.legend.add_face(CircleFace(50, v), column=0)
                    ts.legend.add_face(TextFace(k, fsize=10), column=1)

            # Add title
            ts.title.add_face(TextFace("Colouring tips on " + col, fsize=20), column=0)

    return tree, ts


def get_random_color(pastel_factor=0.5):
    return [
        (x + pastel_factor) / (1.0 + pastel_factor)
        for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]
    ]


def color_distance(c1, c2):
    return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])


def generate_new_color(existing_colors, pastel_factor=0.5):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = get_random_color(pastel_factor=pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color, c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color


def get_color_dict(annot_dict, random_val):
    random.seed(random_val)

    color_dict = {}

    # print (annot_dict)

    # for val in annot_dict.values():
    #     print (val)
    #     print (float(val) is np.nan)
    #     print (np.isnan(val))
    #     print(np.isnan(float(val)))

    unique_vals = list(set(val for val in annot_dict.values()))

    print("uv")
    print("unique_vals")

    colors = []
    for uv in unique_vals:
        color = generate_new_color(colors, pastel_factor=0.9)
        colors.append(color)
        color_dict[uv] = "#%02x%02x%02x" % (
            int(color[0] * 255),
            int(color[1] * 255),
            int(color[2] * 255),
        )
    return color_dict
